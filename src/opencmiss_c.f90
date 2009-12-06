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
  
  PUBLIC CMISSFinaliseC,CMISSInitialiseNumberC,CMISSInitialiseCPtr

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
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT),VALUE, INTENT(IN) :: FieldUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: FileName(FileNameSize)
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisOutputNumberC
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
    TYPE(C_PTR), INTENT(OUT) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: FileName(FileNameSize)
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisOutputPtrC
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
  FUNCTION CMISSInitialiseNumberC(WorldCoordinateSystemUserNumber,WorldRegionUserNumber) BIND(C,NAME="CMISSInitialiseNum")
  
    !Argument variables
    INTEGER(C_INT), INTENT(OUT) :: WorldCoordinateSystemUserNumber
    INTEGER(C_INT), INTENT(OUT) :: WorldRegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSInitialiseNumberC
    !Local variables
  
    CALL CMISSInitialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,CMISSInitialiseNumberC)

    RETURN
    
  END FUNCTION CMISSInitialiseNumberC

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
  FUNCTION CMISSFieldComponentInterpolationGetNumberC(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType) BIND(C,NAME = "CMISSFieldComponentInterpolationGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the interpolation type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the interpolation type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationGetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldComponentInterpolationGetNumber(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType, CMISSFieldComponentInterpolationGetNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationGetNumberC

  !!
  !!==================================================================================================================================
  !!

  !>Returns the interpolation type for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentInterpolationGetPtrC(FieldPtr,VariableType,ComponentNumber,InterpolationType) BIND(C, NAME = &
  & "CMISSFieldComponentInterpolationGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the interpolation type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the interpolation type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentInterpolationGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentInterpolationGetObj(Field, VariableType,ComponentNumber,InterpolationType, &
        & CMISSFieldComponentInterpolationGetPtrC)
      ELSE
        CMISSFieldComponentInterpolationGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentInterpolationGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationGetPtrC

  !!
  !!==================================================================================================================================
  !!

  !>Sets/changes the interpolation type for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentInterpolationSetNumberC(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType) BIND(C, NAME = "CMISSFieldComponentInterpolationSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the interpolation type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldComponentInterpolationSetNumber(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,InterpolationType,&
    & CMISSFieldComponentInterpolationSetNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationSetNumberC


  !!
  !!==================================================================================================================================
  !!

  !>Sets/changes the interpolation type for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentInterpolationSetPtrC(FieldPtr,VariableType,ComponentNumber,InterpolationType) &
   & BIND(C, NAME = "CMISSFieldComponentInterpolationSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the interpolation type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentInterpolationSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentInterpolationSetObj(Field, VariableType, ComponentNumber,InterpolationType, &
        & CMISSFieldComponentInterpolationSetPtrC)
      ELSE
        CMISSFieldComponentInterpolationSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentInterpolationSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationSetPtrC


  !!
  !!==================================================================================================================================
  !!

  !>Returns the character string label for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentLabelGetCNumberC(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,LabelSize,Label) &
  & BIND(C, NAME = "CMISSFieldComponentLabelGetCNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable component character string label to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelGetCNumberC !<Error Code.
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CALL CMISSFieldComponentLabelGetCNumber(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,FLabel, &
    & CMISSFieldComponentLabelGetCNumberC)
    CALL CMISSF2CString(Flabel,Label)

    RETURN

  END FUNCTION CMISSFieldComponentLabelGetCNumberC

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentLabelGetCPtrC(FieldPtr,VariableType,ComponentNumber,LabelSize,Label) BIND(C, NAME = &
  & "CMISSFieldComponentLabelGetCObj")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable component character string label to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelGetCPtrC !<Error Code.
    !Local variable
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentLabelGetCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentLabelGetCObj(Field,VariableType,ComponentNumber,FLabel,CMISSFieldComponentLabelGetCPtrC)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldComponentLabelGetCPtrC=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentLabelGetCPtrC=CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentLabelGetCPtrC

  !
  !================================================================================================================================
  !
  !>Sets/changes the character string label for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentLabelSetCNumberC(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,LabelSize,Label)&
  & BIND(C,NAME= "CMISSFieldComponentLabelSetCNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable component character string label to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelSetCNumberC !<Error Code.
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CALL CMISSC2FString(Label,Flabel)
    CALL CMISSFieldComponentLabelSetCNumber(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,FLabel, &
    & CMISSFieldComponentLabelSetCNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentLabelSetCNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentLabelSetCPtrC(FieldPtr,VariableType,ComponentNumber,LabelSize,Label) BIND(C, NAME= &
  & "CMISSFieldComponentLabelSetCObj")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label Size
    CHARACTER(LEN = 1, KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable component character string label to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelSetCPtrC !<Error Code.
    !Local variable
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentLabelSetCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldComponentLabelSetCObj(Field,VariableType,ComponentNumber,FLabel,CMISSFieldComponentLabelSetCPtrC)
      ELSE
        CMISSFieldComponentLabelSetCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentLabelSetCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSFieldComponentLabelSetCPtrC

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentMeshComponentGetNumberC(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & MeshComponent)

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the mesh component number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the mesh component number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the mesh component number for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the mesh component number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponent !<The mesh component number to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldComponentMeshComponentGetNumber(RegionUserNumber, FieldUserNumber, VariableType, ComponentNumber, MeshComponent,&
    & CMISSFieldComponentMeshComponentGetNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentMeshComponentGetPtrC(FieldPtr,VariableType,ComponentNumber,MeshComponent) BIND(C, &
  & NAME = "CMISSFieldComponentMeshComponentGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the mesh component number for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the mesh component number for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the mesh component number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponent !<The mesh component number to get, for C.
    !Function Variables
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentMeshComponentGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentMeshComponentGetObj(Field, VariableType, ComponentNumber, MeshComponent, &
        & CMISSFieldComponentMeshComponentGetPtrC)
      ELSE
        CMISSFieldComponentMeshComponentGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentMeshComponentGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentGetPtrC

  !
  !================================================================================================================================
  !
  !>Sets/changes the mesh component number for a field variable component for a field identified by a user number.
  FUNCTION CMISSFieldComponentMeshComponentSetNumberC(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & MeshComponent) BIND(C, NAME = "CMISSFieldComponentMeshComponentSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the mesh component number to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the mesh component number to, for C.
    INTEGER(C_INT), INTENT(IN) :: MeshComponent !<The mesh component number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldComponentMeshComponentSetNumber(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,MeshComponent, &
    & CMISSFieldComponentMeshComponentSetNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentMeshComponentSetPtrC(FieldPtr,VariableType,ComponentNumber,MeshComponent) BIND(C, &
  & NAME = "CMISSFieldComponentMeshComponentSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the mesh component number to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the mesh component number to. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponent !<The mesh component number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentMeshComponentSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentMeshComponentSetObj(Field,VariableType,ComponentNumber,MeshComponent, &
        & CMISSFieldComponentMeshComponentSetPtrC)
      ELSE
        CMISSFieldComponentMeshComponentSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentMeshComponentSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentSetPtrC
  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to an integer constant value for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentValuesInitialiseIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseIntgNumber")

        !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseIntgNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldComponentValuesInitialiseIntgNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, &
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseIntgNumberC

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to an integer constant value for a field identified by an object for C.
  FUNCTION CMISSFieldComponentValuesInitialiseIntgPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME ="CMISSFieldComponentValuesInitialiseIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseIntgObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseIntgPtrC)
      ELSE
        CMISSFieldComponentValuesInitialiseIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseIntgPtrC

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a single precision constant value for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentValuesInitialiseSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for , for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseSPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseSPNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, &
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseSPNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseSPNumberC

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a single precision constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseSPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldComponentValuesInitialiseSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseSPPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseSPPtrC)
      ELSE
        CMISSFieldComponentValuesInitialiseSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseSPPtrC

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a double precision constant value for a field identified by a user number, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseDPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseDPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseDPNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType,&
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseDPNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseDPNumberC

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a double precision constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseDPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldComponentValuesInitialiseDP")

      !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseDPPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseDPPtrC)
      ELSE
        CMISSFieldComponentValuesInitialiseDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseDPPtrC


  !>Initialises the values of parameter set of a field variable component to a logical constant value for a field identified by a user number, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseLNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, ComponentNumber,&
    & Value, CMISSFieldComponentValuesInitialiseLNumberC)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseLNumberC

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a logical constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseLPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldComponentValuesInitialiseLObj")

    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseLPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseLObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseLPtrC)
      ELSE
        CMISSFieldComponentValuesInitialiseLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseLPtrC

  !
  !================================================================================================================================
  !

  !>Returns the data type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDataTypeGetNumberC(RegionUserNumber,FieldUserNumber,VariableType,DataType) BIND(C,NAME= &
  & "CMISSFieldDataTypeGetNumber")
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the data type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DataType !<The field variable data type to get, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeGetNumberC !<Error Code.
    !Local variables

    CALL  CMISSFieldDataTypeGetNumber(RegionUserNumber,FieldUserNumber,VariableType,DataType,CMISSFieldDataTypeGetNumberC)

    RETURN

  END FUNCTION  CMISSFieldDataTypeGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the data type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDataTypeGetPtrC(FieldPtr,VariableType,DataType) BIND(C, NAME = "CMISSFieldDataTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<The field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the data type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DataType !<The field variable data type to get, for C. \see OPENCMISS_FieldDataTypes
    !Function Variables
    INTEGER(C_INT) :: CMISSFieldDataTypeGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDataTypeGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDataTypeGetObj(Field, VariableType, DataType, CMISSFieldDataTypeGetPtrC)
      ELSE
        CMISSFieldDataTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDataTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDataTypeGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDataTypeSetNumberC(RegionUserNumber,FieldUserNumber,VariableType,DataType) BIND(C, NAME = &
  & "CMISSFieldDataTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the data type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DataType !<The field variable data type to set, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldDataTypeSetNumber(RegionUserNumber, FieldUserNumber,VariableType, DataType,CMISSFieldDataTypeSetNumberC)

    RETURN

  END FUNCTION CMISSFieldDataTypeSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDataTypeSetPtrC(FieldPtr,VariableType,DataType) BIND(C, NAME = "CMISSFieldDataTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the data type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DataType !<The field variable data type to set, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDataTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDataTypeSetObj(Field, VariableType, DataType,CMISSFieldDataTypeSetPtrC)
      ELSE
        CMISSFieldDataTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDataTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDataTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the DOF order type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDOFOrderTypeGetNumberC(RegionUserNumber,FieldUserNumber,VariableType,DOFOrderType) BIND( &
  & C, NAME = "CMISSFieldDOFOrderTypeGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the DOF order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the DOF order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the DOF order type for, for C.\see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DOFOrderType !<The field variable DOF Order type to get, for C.  \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeGetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldDOFOrderTypeGetNumber(RegionUserNumber, FieldUserNumber, VariableType, DOFOrderType, &
    & CMISSFieldDOFOrderTypeGetNumberC)

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the DOF Order type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDOFOrderTypeGetPtrC(FieldPtr,VariableType,DOFOrderType) BIND(C, NAME = "CMISSFieldDOFOrderTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer of the field to get the DOF Order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the DOF Order type for, for C \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DOFOrderType !<The field variable DOF Order type to get, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDOFOrderTypeGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDOFOrderTypeGetObj(Field, VariableType, DOFOrderType, CMISSFieldDOFOrderTypeGetPtrC)
      ELSE
        CMISSFieldDOFOrderTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDOFOrderTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the DOF order type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldDOFOrderTypeSetNumberC(RegionUserNumber,FieldUserNumber,VariableType,DOFOrderType) BIND(C, &
  & NAME = "CMISSFieldDOFOrderTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the DOF Order type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DOFOrderType !<The field variable DOF Order type to set, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeSetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldDOFOrderTypeSetNumber(RegionUserNumber, FieldUserNumber, VariableType, DOFOrderType, &
    & CMISSFieldDOFOrderTypeSetNumberC)

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the DOF Order type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldDOFOrderTypeSetPtrC(FieldPtr,VariableType,DOFOrderType) BIND(C, NAME = "CMISSFieldDOFOrderTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the DOF Order type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DOFOrderType !<The field variable DOF Order type to set, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDOFOrderTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDOFOrderTypeSetObj(Field, VariableType, DOFOrderType, CMISSFieldDOFOrderTypeSetPtrC)
      ELSE
        CMISSFieldDOFOrderTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDOFOrderTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field identified by a user number for C.
  FUNCTION CMISSFieldCreateFinishNumberC(RegionUserNumber,FieldUserNumber) BIND(C, NAME = "CMISSFieldCreateFinishNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to finish the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to finish the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateFinishNumberC !<Error Code.

    CALL CMISSFieldCreateFinishNumber(RegionUserNumber,FieldUserNumber,CMISSFieldCreateFinishNumberC)

    RETURN

  END FUNCTION CMISSFieldCreateFinishNumberC

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field identified by an object for C.
  FUNCTION CMISSFieldCreateFinishPtrC(FieldPtr) BIND(C, NAME = "CMISSFieldCreateFinishObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finish the creation of.
    !Function variable
    INTEGER(C_INT) ::CMISSFieldCreateFinishPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldCreateFinishPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldCreateFinishObj(Field,CMISSFieldCreateFinishPtrC)
      ELSE
        CMISSFieldCreateFinishPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldCreateFinishPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldCreateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Starts the creation of a field identified by a user number for C.
  FUNCTION CMISSFieldCreateStartNumberC(FieldUserNumber,RegionUserNumber) BIND(C, NAME = "CMISFieldCreateStartNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber!<The user number of the region containing the field to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the field to start the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateStartNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldCreateStartNumber(FieldUserNumber, RegionUserNumber, CMISSFieldCreateStartNumberC)

    RETURN

  END FUNCTION CMISSFieldCreateStartNumberC


  !
  !================================================================================================================================
  !

  !>Starts the creation of a field identified by an object for C.
  FUNCTION CMISSFieldCreateStartPtrC(FieldUserNumber,RegionPtr,FieldPtr) BIND(C, NAME ="CMISSFieldCreateStartObj")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber  !<The user number of the field to start the creation of, for C.
    TYPE(C_PTR), INTENT(IN) :: RegionPtr !<C pointer to the region to create the field on.
    TYPE(C_PTR), INTENT(IN) :: FieldPtr !<C pointer to the created field.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateStartPtrC !<Error Code.
    !Local variable
    TYPE(CMISSRegionType), POINTER :: Region
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldCreateStartPtrC = CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr, Region)
      IF(ASSOCIATED(Region)) THEN
        IF(C_ASSOCIATED(FieldPtr)) THEN
          CALL C_F_POINTER(FieldPtr, Field)
          IF(ASSOCIATED(Field)) THEN
            CALL CMISSFieldCreateStartObj(FieldUserNumber, Region, Field, CMISSFieldCreateStartPtrC)
          ELSE
            CMISSFieldCreateStartPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldCreateStartPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldCreateStartPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldCreateStartPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldCreateStartPtrC


  !
  !================================================================================================================================
  !

  !>Returns the dependent type for a field identified by a user number for C.
  FUNCTION CMISSFieldDependentTypeGetNumberC(RegionUserNumber,FieldUserNumber,DependentType) BIND(C, NAME = &
  & "CMISSFieldDependentTypeGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to get the dependent type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the dependent type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DependentType !<The field dependent type to get, for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldDependentTypeGetNumber(RegionUserNumber, FieldUserNumber, DependentType, CMISSFieldDependentTypeGetNumberC)

    RETURN

  END FUNCTION CMISSFieldDependentTypeGetNumberC


  !
  !================================================================================================================================
  !

  !>Returns the dependent type for a field identified by an object for C.
  FUNCTION CMISSFieldDependentTypeGetPtrC(FieldPtr,DependentType) BIND(C, NAME = "CMISSFieldDependentTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the dependent type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DependentType !<The field dependent type for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeGetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDependentTypeGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF (ASSOCIATED(Field)) THEN
        CALL CMISSFieldDependentTypeGetObj(Field, DependentType, CMISSFieldDependentTypeGetPtrC)
      ELSE
        CMISSFieldDependentTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDependentTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDependentTypeGetPtrC


  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field identified by a user number for C.
  FUNCTION CMISSFieldDependentTypeSetNumberC(RegionUserNumber,FieldUserNumber,DependentType) BIND(C, NAME  = &
  & "CMISSFieldDependentTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DependentType !<The field dependent type for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldDependentTypeSetNumber(RegionUserNumber,FieldUserNumber,DependentType,CMISSFieldDependentTypeSetNumberC)

    RETURN

  END FUNCTION CMISSFieldDependentTypeSetNumberC


  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field identified by an object for C.
  FUNCTION CMISSFieldDependentTypeSetPtrC(FieldPtr,DependentType) BIND(C, NAME = "CMISSFieldDependentTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DependentType !<The field dependent type, for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDependentTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDependentTypeSetObj(Field, DependentType, CMISSFieldDependentTypeSetPtrC)
      ELSE
        CMISSFieldDependentTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDependentTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDependentTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Destroys a field identified by a user number for C.
  FUNCTION CMISSFieldDestroyNumberC(RegionUserNumber,FieldUserNumber) BIND(C, NAME = "CMISSFieldDestroyNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to destroy for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDestroyNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldDestroyNumber(RegionUserNumber, FieldUserNumber, CMISSFieldDestroyNumberC)

    RETURN

  END FUNCTION CMISSFieldDestroyNumberC

  !
  !================================================================================================================================
  !

  !>Destroys a field identified by an object for C.
  FUNCTION CMISSFieldDestroyPtrC(FieldPtr) BIND(C,NAME = "CMISSFieldDestroyObj")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDestroyPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDestroyPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF (ASSOCIATED(Field)) THEN
        CALL CMISSFieldDestroyObj(Field, CMISSFieldDestroyPtrC)
      ELSE
        CMISSFieldDestroyPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDestroyPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDestroyPtrC

  !
  !================================================================================================================================
  !

  !>Returns the dimension for a field identified by a user number for C.
  FUNCTION CMISSFieldDimensionGetNumberC(RegionUserNumber,FieldUserNumber,VariableType,DIMENSION) BIND(C, NAME = &
  & "CMISSFieldDimensionGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the dimension for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionGetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldDimensionGetNumber(RegionUserNumber, FieldUserNumber, VariableType, DIMENSION, CMISSFieldDimensionGetNumberC)

    RETURN

  END FUNCTION CMISSFieldDimensionGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the dimension for a field identified by an object for C.
  FUNCTION CMISSFieldDimensionGetPtrC(FieldPtr,VariableType,Dimension) BIND(C, NAME = "CMISSFieldDimensionGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the dimension for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionGetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDimensionGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDimensionGetObj(Field, VariableType, Dimension, CMISSFieldDimensionGetPtrC)
      ELSE
        CMISSFieldDimensionGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDimensionGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDimensionGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension for a field identified by a user number for C.
  FUNCTION CMISSFieldDimensionSetNumberC(RegionUserNumber,FieldUserNumber,VariableType,Dimension) BIND(C, NAME = &
  & "CMISSFieldDimensionSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the dimension to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionSetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldDimensionSetNumber(RegionUserNumber, FieldUserNumber, VariableType, Dimension, CMISSFieldDimensionSetNumberC)

    RETURN

  END FUNCTION CMISSFieldDimensionSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension for a field identified by an object for C.
  FUNCTION CMISSFieldDimensionSetPtrC(FieldPtr,VariableType,Dimension) BIND(C, NAME= "CMISSFieldDimensionSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the dimension to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: Dimension !<The field dimension, for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDimensionSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDimensionSetObj(Field, VariableType, Dimension, CMISSFieldDimensionSetPtrC)
      ELSE
        CMISSFieldDimensionSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDimensionSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDimensionSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the geometric field for a field identified by a user number for C.
  FUNCTION CMISSFieldGeometricFieldGetNumberC(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber) BIND(C, NAME = &
  & "CMISSFieldGeometricFieldGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the geometric field for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the geometric field for, for C.
    INTEGER(C_INT), INTENT(OUT) :: GeometricFieldUserNumber !<The field geometric field user number, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldGeometricFieldGetNumber(RegionUserNumber, FieldUserNumber, GeometricFieldUserNumber, &
    & CMISSFieldGeometricFieldGetNumberC)

    RETURN

  END FUNCTION CMISSFieldGeometricFieldGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the geometric field for a field identified by an object for C.
  FUNCTION CMISSFieldGeometricFieldGetPtrC(FieldPtr,GeometricFieldPtr) BIND(C, NAME = "CMISSFieldGeometricFieldGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the geometric field to, for C.
    TYPE(C_PTR), INTENT(OUT) :: GeometricFieldPtr !<C pointer to the geometric field for the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSFieldType), POINTER :: GeometricField

    CMISSFieldGeometricFieldGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeometricFieldPtr)) THEN
          CALL C_F_POINTER(GeometricFieldPtr, GeometricField)
          IF(ASSOCIATED(GeometricField)) THEN
            CALL CMISSFieldGeometricFieldGetObj(Field, GeometricField, CMISSFieldGeometricFieldGetPtrC)
          ELSE
            CMISSFieldGeometricFieldGetPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldGeometricFieldGetPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldGeometricFieldGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldGeometricFieldGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldGeometricFieldGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field identified by a user number for C.
  FUNCTION CMISSFieldGeometricFieldSetNumberC(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber) BIND(C, &
  & NAME = "CMISSFieldGeometricFieldSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region corresponding to the field to set the geometric field to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !< The user number for the field to set the geometric field to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeometricFieldUserNumber !<The field geometric field user number to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldSetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldGeometricFieldSetNumber(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber,&
    & CMISSFieldGeometricFieldSetNumberC)

    RETURN

  END FUNCTION CMISSFieldGeometricFieldSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field identified by an object for C.
  FUNCTION CMISSFieldGeometricFieldSetPtrC(FieldPtr,GeometricFieldPtr) BIND(C, NAME = "CMISSFieldGeometricFieldSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the geometric field to, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeometricFieldPtr !<C pointer to the geometric field for the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSFieldType), POINTER :: GeometricField

    CMISSFieldGeometricFieldSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeometricFieldPtr)) THEN
          CALL C_F_POINTER(GeometricFieldPtr, GeometricField)
          IF(ASSOCIATED(GeometricField)) THEN
            CALL CMISSFieldGeometricFieldSetObj(Field, GeometricField, CMISSFieldGeometricFieldSetPtrC)
          ELSE
            CMISSFieldGeometricFieldSetPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldGeometricFieldSetPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldGeometricFieldSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
     CMISSFieldGeometricFieldSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldGeometricFieldSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field identified by a user number for C.
  FUNCTION CMISSFieldLabelGetCNumberC(RegionUserNumber,FieldUserNumber,LabelSize,Label) BIND(C, NAME = "CMISSFieldLabelGetCNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label Size
    CHARACTER(LEN = 1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelGetCNumberC !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize - 1) :: FLabel

    CALL CMISSFieldLabelGetCNum(RegionUserNumber,FieldUserNumber,FLabel,CMISSFieldLabelGetCNumberC)
    CALL CMISSF2CString(FLabel, Label)

    RETURN

  END FUNCTION CMISSFieldLabelGetCNumberC

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field identified by an object for C.
  FUNCTION CMISSFieldLabelGetCPtrC(FieldPtr,LabelSize,Label) BIND(C, NAME = "CMISSFieldLabelGetC")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label size
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelGetCPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize -1) :: FLabel

    CMISSFieldLabelGetCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldLabelGetC(Field, FLabel, CMISSFieldLabelGetCPtrC)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldLabelGetCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldLabelGetCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldLabelGetCPtrC

  !!
  !================================================================================================================================
  !!

  !>Sets/changes the character string label for a field identified by a user number for C.
  FUNCTION CMISSFieldLabelSetCNumberC(RegionUserNumber,FieldUserNumber,LabelSize, Label) BIND(C, NAME = "CMISSFieldLabelSetCNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size.
    CHARACTER(LEN = 1, KIND= C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelSetCNumberC !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSC2FString (Label, FLabel)
    CALL CMISSFieldLabelSetCNumber(RegionUserNumber, FieldUserNumber, FLabel,CMISSFieldLabelSetCNumberC)

    RETURN

  END FUNCTION CMISSFieldLabelSetCNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field identified by an object for C.
  FUNCTION CMISSFieldLabelSetCPtrC(FieldPtr,LabelSize, Label) BIND(C, NAME = "CMISSFieldLabelSetCObj")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelSetCPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CMISSFieldLabelSetCPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldLabelSetCObj(Field,FLabel,CMISSFieldLabelSetCPtrC)
      ELSE
        CMISSFieldLabelSetCPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldLabelSetCPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldLabelSetCPtrC


  !
  !================================================================================================================================
  !

  !>Returns the mesh decomposition for a field identified by a user number for C.

  FUNCTION CMISSFieldMeshDecompositionGetNumberC(RegionUserNumber,FieldUserNumber,DecompositionUserNumber) &
    & BIND(C, NAME = "CMISSFieldMeshDecompositionGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The field mesh decomposition user number for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldMeshDecompositionGetNumber(RegionUserNumber,FieldUserNumber,DecompositionUserNumber,&
    & CMISSFieldMeshDecompositionGetNumberC)

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the mesh decomposition for a field identified by an object for C.
  FUNCTION CMISSFieldMeshDecompositionGetPtrC(FieldPtr,MeshDecompositionPtr) BIND(C, NAME = "CMISSFieldMeshDecompositionGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the mesh decomposition for, for C.
    TYPE(C_PTR), INTENT(OUT) :: MeshDecompositionPtr !<C pointer to the field mesh decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSDecompositionType), POINTER :: MeshDecomposition

    CMISSFieldMeshDecompositionGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(MeshDecompositionPtr)) THEN
          CALL C_F_POINTER(MeshDecompositionPtr, MeshDecomposition)
          IF(ASSOCIATED(MeshDecomposition)) THEN
            CALL CMISSFieldMeshDecompositionGetObj(Field, MeshDecomposition, CMISSFieldMeshDecompositionGetPtrC)
          ELSE
            CMISSFieldMeshDecompositionGetPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldMeshDecompositionGetPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldMeshDecompositionGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldMeshDecompositionGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field identified by a user number for C.
  FUNCTION CMISSFieldMeshDecompositionSetNumberC(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber) &
  & BIND(C, NAME = "CMISSFieldMeshDecompositionSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number for the mesh to set the mesh decomposititon to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The field mesh decomposition user number for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldMeshDecompositionSetNumber(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber,&
    & CMISSFieldMeshDecompositionSetNumberC)

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field identified by an object for C.
  FUNCTION CMISSFieldMeshDecompositionSetPtrC(FieldPtr,MeshDecompositionPtr) BIND(C,NAME = "CMISSFieldMeshDecompositionSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the mesh decomposition to, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshDecompositionPtr !<C pointer to the mesh decomposition for the field to set.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSDecompositionType), POINTER :: MeshDecomposition

    CMISSFieldMeshDecompositionSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(MeshDecompositionPtr)) THEN
          CALL C_F_POINTER(MeshDecompositionPtr, MeshDecomposition)
          IF (ASSOCIATED(MeshDecomposition)) THEN
            CALL CMISSFieldMeshDecompositionSetObj(Field, MeshDecomposition, CMISSFieldMeshDecompositionSetPtrC)
          ELSE
            CMISSFieldMeshDecompositionSetPtrC = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldMeshDecompositionSetPtrC = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldMeshDecompositionSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldMeshDecompositionSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the number of componenets for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfComponentsGetNumberC(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents) &
  & BIND(C, NAME = "CMISSFieldNumberOfComponentsGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsGetNumberC
    !Local variables

    CALL CMISSFieldnumberOfComponentsGetNumber(RegionUserNumber, FieldUserNumber, VariableType, NumberOfComponents, &
    & CMISSFieldNumberOfComponentsGetNumberC)

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the number of components for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfComponentsGetPtrC(FieldPtr,VariableType,NumberOfComponents) BIND(C, NAME = &
  & "CMISSFieldNumberOfComponentsGetObj")

    !Argument variables
    TYPE(C_PTR),VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the number of components for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the field variable for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfComponentsGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfComponentsGetObj(Field, VariableType, NumberOfComponents)
      ELSE
        CMISSFieldNumberOfComponentsGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfComponentsGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsGetPtrC


  !
  !================================================================================================================================
  !

  !>Sets/changes the number of componenets for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfComponentsSetNumberC(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents) &
  & BIND(C, NAME = "CMISSFieldNumberOfComponentsSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the number of components to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the field variable for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldNumberOfComponentsSetNumber(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents, &
    & CMISSFieldNumberofComponentsSetNumberC)

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfComponentsSetPtrC(FieldPtr,VariableType,NumberOfComponents) BIND(C, NAME = &
  & "CMISSFieldNumberOfComponentsSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the number of components for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the field variables for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfComponentsSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfComponentsSetObj(Field, VariableType, NumberOfComponents, CMISSFieldNumberOfComponentsSetPtrC)
      ELSE
        CMISSFieldNumberOfComponentsSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfComponentsSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the number of variables for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfVariablesGetNumberC(RegionUserNumber,FieldUserNumber,NumberOfVariables) BIND(C, NAME = &
  & "CMISSFieldNumberOfVariablesGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the number of variables for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the number of variables for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfVariables !<The number of variables in the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldNumberOfVariablesGetNumber(RegionUserNumber, FieldUserNumber, NumberOfVariables, &
    & CMISSFieldNumberOfVariablesGetNumberC)

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the number of variables for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfVariablesGetPtrC(FieldPtr,NumberOfVariables) BIND(C, NAME = "CMISSFieldNumberOfVariablesGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the number of variables for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfVariables !<The number of variables in the field for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfVariablesGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfVariablesGetObj(Field, NumberOfVariables, CMISSFieldNumberOfVariablesGetPtrC)
      ELSE
        CMISSFieldNumberOfVariablesGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfVariablesGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfVariablesSetNumberC(RegionUserNumber,FieldUserNumber,NumberOfVariables) BIND(C, &
  & NAME = "CMISSFieldNumberOfVariablesSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the number of variables to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the number of variables to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfVariables !<The number of variables set to the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesSetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldNumberofVariablesSetNumber(RegionUserNumber, FieldUserNumber, NumberOfVariables, &
    & CMISSFieldNumberOfVariablesSetNumberC)

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfVariablesSetPtrC(FieldPtr,NumberOfVariables) BIND(C, NAME = "CMISSFieldNumberOfVariablesSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfVariables
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesSetPtrC
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfVariablesSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfVariablesSetObj(Field, NumberOfVariables, CMISSFieldNumberOfVariablesSetPtrC)
      ELSE
        CMISSFieldNumberOfVariablesSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfVariablesSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesSetPtrC


  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME= "CMISSFieldParameterSetAddConstantIntgNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantIntgNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantIntgNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value, CMISSFieldParameterSetAddConstantIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantIntgNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddConstantIntgPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME= "CMISSFieldParameterSetAddConstantIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantIntgPtrC
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantIntgPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantIntgObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantIntgPtrC)
      ELSE
        CMISSFieldParameterSetAddConstantIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantIntgPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.

  FUNCTION CMISSFieldParameterSetAddConstantSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantSPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantSPNumber (RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetAddConstantSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantSPNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddConstantSPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetAddConstantSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantSPObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantSPPtrC)
      ELSE
        CMISSFieldParameterSetAddConstantSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantSPPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantDPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantDPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantDPNumber (RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetAddConstantDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantDPNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddConstantDPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetAddConstantDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantDPObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantDPPtrC)
      ELSE
        CMISSFieldParameterSetAddConstantDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantDPPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetAddConstantLNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, ComponentNumber, &
    & Value, CMISSFieldParameterSetAddConstantLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantLNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddConstantLPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetAddConstantLPtrC")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantLPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantLObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantLPtrC)
      ELSE
        CMISSFieldParameterSetAddConstantLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN


  END FUNCTION CMISSFieldParameterSetAddConstantLPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementIntgNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementIntgNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementIntgNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value, CMISSFieldParameterSetAddElementIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementIntgNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an element in the given parameter set for field variable component for a field identified by an object, for C.
  FUNCTION CMISSFieldParameterSetAddElementIntgPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementIntgPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementIntgObj(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementIntgPtrC)
      ELSE
        CMISSFieldParameterSetAddElementIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementIntgPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an element in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddElementSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementSPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementSPNumber (RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementSPNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an element in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddElementSPPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementSPPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementSPObj(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementSPPtrC)
      ELSE
        CMISSFieldParameterSetAddElementSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementSPPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementDPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementDPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementDPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementDPNumberC


    !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an element in the given parameter set for field variable component for a field identified by an object, for C.
  FUNCTION CMISSFieldParameterSetAddElementDPPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementDPPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementDPObj(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementDPPtrC)
      ELSE
        CMISSFieldParameterSetAddElementDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementDPPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementLNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementLNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementLNumberC

    !
  !================================================================================================================================
  !

  !>Adds the given logical value to an element in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddElementLPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementLObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementLPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementLPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementLObj(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
         & CMISSFieldParameterSetAddElementLPtrC)
      ELSE
        CMISSFieldParameterSetAddElementLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementLPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an node in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddNodeIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
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
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeIntgNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeIntgNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber,&
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeIntgNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeIntgPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeIntgPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeIntgObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeIntgPtrC)
      ELSE
        CMISSFieldParameterSetAddNodeIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeIntgPtrC



  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeSPNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeSPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeSPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeSPNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeSPPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value  !<The single precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeSPPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeSPObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeSPPtrC)
      ELSE
        CMISSFieldParameterSetAddNodeSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeSPPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeDPNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeDPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeDPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeDPNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an node in the given parameter set for field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddNodeDPPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value  !<The double precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeDPPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeDPObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeDPPtrC)
      ELSE
        CMISSFieldParameterSetAddNodeDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeDPPtrC

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeLNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeLNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeLNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddNodeLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeLNumberC

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeLPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber,ComponentNumber, &
    & Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeLObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeLPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeLPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeLObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddNodeLPtrC)
      ELSE
        CMISSFieldParameterSetAddNodeLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeLPtrC

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetCreateNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, &
  & NAME = "CMISSFieldParameterSetCreateNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to create the parameter set on for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to create, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetCreateNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetCreateNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetCreateNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetCreateNumberC

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetCreatePtrC(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetCreateObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to create the parameter set on for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to create, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetCreatePtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetCreatePtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetCreateObj(Field,VariableType,FieldSetType,CMISSFieldParameterSetCreatePtrC)
      ELSE
        CMISSFieldParameterSetCreatePtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetCreatePtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetCreatePtrC

  !
  !================================================================================================================================
  !

  !>Destroys the specified parameter set type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetDestroyNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, NAME= &
  & "CMISSFieldParameterSetDestroyNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to destroy the parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to destroy, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDestroyNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetDestroyNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetDestroyNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetDestroyNumberC

  !
  !================================================================================================================================
  !

  !>Destroys the specified parameter set type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetDestroyPtrC(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetDestroyObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to destroy the parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to destroy, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDestroyPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetDestroyPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetDestroyObj(Field,VariableType,FieldSetType,CMISSFieldParameterSetDestroyPtrC)
      ELSE
        CMISSFieldParameterSetDestroyPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetDestroyPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDestroyPtrC

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local integer data array for a field identified by an user number, for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersPtr) &
  & BIND(C, NAME = "CMISSFieldParameterSetDataGetIntgNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get for C. \see OPENCMISS_FieldParameterSetTypes
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetIntgNumberC
    !Local variables

    CMISSFieldParameterSetDataGetIntgNumberC = CMISSNoError
    CALL CMISSFieldParameterSetDataGetIntgNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType, ParametersPtr, &
    & CMISSFieldParameterSetDataGetIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetIntgNumberC

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local integer data array for a field identified by an object for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetIntgPtrC(FieldPtr,VariableType,FieldSetType,ParametersPtr) BIND(C, NAME = &
  & "CMISSFieldParameterSetDataGetIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetIntgPtrC
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetDataGetIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetDataGetIntgObj(Field, VariableType, FieldSetType, ParametersPtr, &
        & CMISSFieldParameterSetDataGetIntgPtrC)
      ELSE
        CMISSFieldParameterSetDataGetIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetDataGetIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetIntgPtrC

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local single precision data array for a field identified by an user number for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
  &ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataGetSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetSPNumberC
    !Local variables

    CMISSFieldParameterSetDataGetSPNumberC = CMISSNoError
    CALL CMISSFieldParameterSetDataGetSPNumber(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, ParametersPtr, &
    & CMISSFieldParameterSetDataGetSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetSPNumberC




!missing code


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantIntgNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantIntgNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantIntgNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,Value, &
    & CMISSFieldParameterSetGetConstantIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantIntgNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantIntgPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value)

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantIntgPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantIntgObj(Field,VariableType,FieldSetType,ComponentNumber,Value, &
        & CMISSFieldParameterSetGetConstantIntgPtrC)
      ELSE
        CMISSFieldParameterSetGetConstantIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantIntgPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantSPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantSPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
    & Value,CMISSFieldParameterSetGetConstantSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantSPNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantSPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetGetConstantSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantSPPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantSPObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetGetConstantSPPtrC)
      ELSE
        CMISSFieldParameterSetGetConstantSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantSPPtrC


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantDPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantDPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantDPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
    & Value,CMISSFieldParameterSetGetConstantDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantDPNumberC

    !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantDPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetGetConstantDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantDPPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantDPObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetGetConstantDPPtrC)
      ELSE
        CMISSFieldParameterSetGetConstantDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantDPPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantLNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
    & Value,CMISSFieldParameterSetGetConstantLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantLNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantLPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetGetConstantLObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantLPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantLObj(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetGetConstantLPtrC)
      ELSE
        CMISSFieldParameterSetGetConstantLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantLPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
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
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementIntgNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementIntgNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetGetElementIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementIntgNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementIntgPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetGetElementIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementIntgObj(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetGetElementIntgPtrC)
      ELSE
        CMISSFieldParameterSetGetElementIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementIntgPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementSPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementSPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetGetElementSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementSPNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified element of a field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetGetElementSPPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetGetElementSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementSPObj(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetGetElementSPPtrC)
      ELSE
        CMISSFieldParameterSetGetElementSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementSPPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementDPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
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
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementDPNumber !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementDPNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetGetElementDPNumber)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementDPNumber

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementDPPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetGetElementDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementDPObj(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetGetElementDPPtrC)
      ELSE
        CMISSFieldParameterSetGetElementDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementDPPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementLNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    &ComponentNumber,Value,CMISSFieldParameterSetGetElementLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementLNumberC


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementLPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetGetElementLObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementLPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementLObj(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value,  &
        & CMISSFieldParameterSetGetElementLPtrC)
      ELSE
        CMISSFieldParameterSetGetElementLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementLPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeIntgNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeIntgNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeIntgNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
    & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeIntgNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeIntgPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeIntgPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeIntgObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetGetNodeIntgPtrC)
      ELSE
        CMISSFieldParameterSetGetNodeIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeIntgPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeSPNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeSPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeSPNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
    &UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeSPNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeSPPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeSPPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeSPObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetGetNodeSPPtrC)
      ELSE
        CMISSFieldParameterSetGetNodeSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeSPPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeDPNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeDPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeDPNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
    & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeDPNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeDPPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeDPPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeDPObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetGetNodeDPPtrC)
      ELSE
        CMISSFieldParameterSetGetNodeDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeDPPtrC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeLNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeLNumber(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
    & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeLNumberC

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeLPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeLObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeLPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeLObj(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetGetNodeLPtrC)
      ELSE
        CMISSFieldParameterSetGetNodeLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeLPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantIntgNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantIntgNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantIntgNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
    & ComponentNumber,Value,CMISSFieldParameterSetUpdateConstantIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantIntgNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetUpdateConstantIntgPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantIntgObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantIntgObj(Field, VariableType,FieldSetType,ComponentNumber,Value,&
        & CMISSFieldParameterSetUpdateConstantIntgPtrC)
      ELSE
        CMISSFieldParameterSetUpdateConstantIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantIntgPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantSPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantSPNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetUpdateConstantSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantSPNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantSPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetUpdateConstantSPObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantSPObj(Field, VariableType,FieldSetType,ComponentNumber,Value, &
        & CMISSFieldParameterSetUpdateConstantSPPtrC)
      ELSE
        CMISSFieldParameterSetUpdateConstantSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantSPPtrC
  
  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantDPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantDPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantDPNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
    & ComponentNumber,Value,CMISSFieldParameterSetUpdateConstantDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantDPNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantDPPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetUpdateConstantDPObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantDPObj(Field, VariableType,FieldSetType,ComponentNumber,Value,&
        & CMISSFieldParameterSetUpdateConstantDPPtrC)
      ELSE
        CMISSFieldParameterSetUpdateConstantDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantDPPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantLNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
    & Value,CMISSFieldParameterSetUpdateConstantLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantLNumberC


  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantLPtrC(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetUpdateConstantLObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantLPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantLObj(Field, VariableType,FieldSetType,ComponentNumber,Value,&
        & CMISSFieldParameterSetUpdateConstantLPtrC)
      ELSE
        CMISSFieldParameterSetUpdateConstantLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantLPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementIntgNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    INTEGER(C_INT), INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementIntgNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementIntgNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
    & UserElementNumber,ComponentNumber,Value,CMISSFieldParameterSetUpdateElementIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementIntgNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementIntgPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementIntgObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    INTEGER(C_INT), INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementIntgObj(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value,&
        & CMISSFieldParameterSetUpdateElementIntgPtrC)
      ELSE
        CMISSFieldParameterSetUpdateElementIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementIntgPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementSPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementSPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementSPNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementSPNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementSPPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementSPObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementSPObj(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value, &
        & CMISSFieldParameterSetUpdateElementSPPtrC)
      ELSE
        CMISSFieldParameterSetUpdateElementSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementSPPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementDPNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementDPNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementDPNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementDPNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementDPPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementDPObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementDPObj(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value,&
        & CMISSFieldParameterSetUpdateElementDPPtrC)
      ELSE
        CMISSFieldParameterSetUpdateElementDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementDPPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementLNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementLNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementLNumber(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementLNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementLPtrC(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementLObj")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementLPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementLObj(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value, &
        & CMISSFieldParameterSetUpdateElementLPtrC)
      ELSE
        CMISSFieldParameterSetUpdateElementLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementLPtrC

  !
  !================================================================================================================================
  !

  !>Finishes the parameter set update for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateFinishNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetUpdateFinishNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to finish the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to finish the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to finish the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to finish the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateFinishNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateFinishNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetUpdateFinishNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateFinishNumberC

  !
  !================================================================================================================================
  !

  !>Finishes the parameter set update for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateFinishPtrC(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = &
  & "CMISSFieldParameterSetUpdateFinishObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to finish the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to finish the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateFinishPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateFinishPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateFinishObj(Field,VariableType,FieldSetType,CMISSFieldParameterSetUpdateFinishPtrC)
      ELSE
        CMISSFieldParameterSetUpdateFinishPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateFinishPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateFinishPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeIntgNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeIntgNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeIntgNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeIntgNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
    & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeIntgNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeIntgNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeIntgPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeIntgObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeIntgPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeIntgPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeIntgObj(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber,Value,CMISSFieldParameterSetUpdateNodeIntgPtrC)
      ELSE
        CMISSFieldParameterSetUpdateNodeIntgPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeIntgPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeIntgPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeSPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeSPNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeSPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeSPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
    & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeSPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeSPNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the node and derivative of the field variable component for a field identified by an object, for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeSPPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeSPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeSPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeSPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeSPObj(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeSPPtrC)
      ELSE
        CMISSFieldParameterSetUpdateNodeSPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeSPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeSPPtrC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeDPNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeDPNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeDPNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeDPNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
    & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeDPNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeDPNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeDPPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeDPObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeDPPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeDPPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeDPObj(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeDPPtrC)
      ELSE
        CMISSFieldParameterSetUpdateNodeDPPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeDPPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeDPPtrC

  !
  !================================================================================================================================
  !,

  !>Updates the given parameter set with the given logical value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeLNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeLNumber")

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
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeLNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeLNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber,&
    & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeLNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeLNumberC

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeLPtrC(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeLObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeLPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeLPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeLObj(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeLPtrC)
      ELSE
        CMISSFieldParameterSetUpdateNodeLPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeLPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeLPtrC

  !
  !================================================================================================================================
  !

  !>Starts the parameter set update for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateStartNumberC(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, &
  & NAME = "CMISSFieldParameterSetUpdateStartNumber")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber!<The user number of the region containing the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to start the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to start the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateStartNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateStartNumber(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,&
    & CMISSFieldParameterSetUpdateStartNumberC)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateStartNumberC

  !
  !================================================================================================================================
  !

  !>Starts the parameter set update for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateStartPtrC(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = &
  & "CMISSFieldParameterSetUpdateStartObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to start the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to start the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateStartPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateStartPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateStartObj(Field,VariableType,FieldSetType,CMISSFieldParameterSetUpdateStartPtrC)
      ELSE
        CMISSFieldParameterSetUpdateStartPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateStartPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateStartPtrC

  !
  !================================================================================================================================
  !

  !>Returns the scaling type for a field identified by a user number for C.
  FUNCTION CMISSFieldScalingTypeGetNumberC(RegionUserNumber,FieldUserNumber,ScalingType) BIND(C, NAME = &
  & "CMISSFieldScalingTypeGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the scaling type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the scaling type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ScalingType !<The field scaling type to get, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeGetNumberC !<Error Code.
    !Local variables

    CALL CMISSFieldScalingTypeGetNumber(RegionUserNumber,FieldUserNumber,ScalingType,CMISSFieldScalingTypeGetNumberC)

    RETURN

  END FUNCTION CMISSFieldScalingTypeGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the scaling type for a field identified by an object for C.
  FUNCTION CMISSFieldScalingTypeGetPtrC(FieldPtr,ScalingType) BIND(C, NAME = "CMISSFieldScalingTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the scaling type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ScalingType !<The field scaling type to get, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeGetPtrC
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldScalingTypeGetPtrC =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldScalingTypeGetObj(Field, ScalingType, CMISSFieldScalingTypeGetPtrC)
      ELSE
        CMISSFieldScalingTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldScalingTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldScalingTypeGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field identified by a user number for C.
  FUNCTION CMISSFieldScalingTypeSetNumberC(RegionUserNumber,FieldUserNumber,ScalingType) BIND(C, NAME = &
  & "CMISSFieldScalingTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ScalingType !<The field scaling type to set, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeSetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldScalingTypeSetNumber(RegionUserNumber,FieldUserNumber,ScalingType,CMISSFieldScalingTypeSetNumberC)

    RETURN

  END FUNCTION CMISSFieldScalingTypeSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field identified by an object for C.
  FUNCTION CMISSFieldScalingTypeSetPtrC(FieldPtr,ScalingType) BIND(C, NAME = "CMISSFieldScalingTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ScalingType !<The field scaling type to set, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldScalingTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldScalingTypeSetObj(Field, ScalingType, CMISSFieldScalingTypeSetPtrC)
      ELSE
        CMISSFieldScalingTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldScalingTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldScalingTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the field type for a field identified by a user number for C.
  FUNCTION CMISSFieldTypeGetNumberC(RegionUserNumber,FieldUserNumber,FieldType) BIND(C, NAME = "CMISSFieldTypeGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the field type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the field type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to get, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeGetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldTypeGetNumber(RegionUserNumber,FieldUserNumber,FieldType,CMISSFieldTypeGetNumberC)

    RETURN

  END FUNCTION CMISSFieldTypeGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the type for a field identified by an object for C.
  FUNCTION CMISSFieldTypeGetPtrC(FieldPtr,FieldType) BIND(C, NAME = "CMISSFieldTypeGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to get, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeGetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldTypeGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldTypeGetObj(Field, FieldType, CMISSFieldTypeGetPtrC)
      ELSE
        CMISSFieldTypeGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the field type for a field identified by a user number for C.
  FUNCTION CMISSFieldTypeSetNumberC(RegionUserNumber,FieldUserNumber,FieldType) BIND(C, NAME = "CMISSFieldTypeSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to set, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeSetNumberC !<Error Code.
    !Local variable

    CALL CMISSFieldTypeSetNumber(RegionUserNumber,FieldUserNumber,FieldType,CMISSFieldTypeSetNumberC)

    RETURN

  END FUNCTION CMISSFieldTypeSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a field identified by an object for C.
  FUNCTION CMISSFieldTypeSetPtrC(FieldPtr,FieldType) BIND(C, NAME = "CMISSFieldTypeSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to set, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeSetPtrC !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldTypeSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldTypeSetObj(Field, FieldType, CMISSFieldTypeSetPtrC)
      ELSE
        CMISSFieldTypeSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable for a field identified by a user number.
  FUNCTION CMISSFieldVariableLabelGetCNumberC(RegionUserNumber,FieldUserNumber,VariableType,LabelSize,Label) BIND(C, &
  & NAME = "CMISSFieldVariableLabelGetCNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelGetCNumberC !<Error Code.
    !Local variable
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSFieldVariableLabelGetCNumber(RegionUserNumber,FieldUserNumber,VariableType,FLabel,CMISSFieldVariableLabelGetCNumberC)
    CALL CMISSF2CString(FLabel,Label)

    RETURN

  END FUNCTION CMISSFieldVariableLabelGetCNumberC

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
    !Local variable

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

  !>Returns the extent for a generated mesh identified by a user number.
  FUNCTION CMISSGeneratedMeshExtentGetNumberC(GeneratedMeshUserNumber,Extent) BIND(C, NAME = "CMISSGeneratedMeshExtentGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to get the extent for, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Extent(:) !<Extent(i). The extent for the i'th dimension of the generated mesh to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshExtentGetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshExtentGetNumber(GeneratedMeshUserNumber, Extent,CMISSGeneratedMeshExtentGetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshExtentGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the extent for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshExtentGetPtrC(GeneratedMeshPtr,Extent) BIND(C, NAME = "CMISSGeneratedMeshExtentGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to get the extent for.
    REAL(C_DOUBLE), INTENT(OUT) :: Extent(:) !<Extent(i). The extent for the i'th dimension of the generated mesh to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshExtentGetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshExtentGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshExtentGetObj(GeneratedMesh, Extent, CMISSGeneratedMeshExtentGetPtrC)
      ELSE
        CMISSGeneratedMeshExtentGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshExtentGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshExtentGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the extent for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshExtentSetNumberC(GeneratedMeshUserNumber,Extent) BIND(C, NAME = "CMISSGeneratedMeshExtentSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the extent to, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Extent(:) !<Extent(i). On return, the extent for the i'th dimension of the generated mesh to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshExtentSetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshExtentSetNumber(GeneratedMeshUserNumber, Extent,CMISSGeneratedMeshExtentSetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshExtentSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the extent for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshExtentSetPtrC(GeneratedMeshPtr,Extent) BIND(C, NAME = "CMISSGeneratedMeshExtentSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to set the extent to.
    REAL(C_DOUBLE), INTENT(OUT) :: Extent(:) !<Extent(i). The extent for the i'th dimension of the generated mesh to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshExtentSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshExtentSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshExtentSetObj(GeneratedMesh, Extent, CMISSGeneratedMeshExtentSetPtrC)
      ELSE
        CMISSGeneratedMeshExtentSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshExtentSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshExtentSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the number of elements for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshNumberOfElementsGetNumberC(GeneratedMeshUserNumber,NumberOfElements) BIND(C, NAME = &
  & "CMISSGeneratedMeshNumberOfElementsGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the number of elements for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements(:) !<NumberOfElements(i). On return, the number of elements in the i'th dimension of the generated mesh, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshNumberOfElementsGetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshNumberOfElementsGetNumber(GeneratedMeshUserNumber, NumberOfElements, &
    & CMISSGeneratedMeshNumberOfElementsGetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshNumberOfElementsGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the number of elements for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshNumberOfElementsGetPtrC(GeneratedMeshPtr,NumberOfElements) BIND(C, NAME = &
  & "CMISSGeneratedMeshNumberOfElementsGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to get the number of elements for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements(:) !<NumberOfElements(i). On return, the number of elements in the i'th dimension of the generated mesh, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshNumberOfElementsGetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshNumberOfElementsGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshNumberOfElementsGetObj(GeneratedMesh, NumberOfElements, CMISSGeneratedMeshNumberOfElementsGetPtrC)
      ELSE
        CMISSGeneratedMeshNumberOfElementsGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshNumberOfElementsGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshNumberOfElementsGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshNumberOfElementsSetNumberC(GeneratedMeshUserNumber,NumberOfElements) BIND(C, NAME = &
  & "CMISSGeneratedMeshNumberOfElementsSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the number of elements to, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements(:) !<NumberOfElements(i). On return, the number of elements in the i'th dimension of the generated mesh to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshNumberOfElementsSetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshNumberOfElementsSetNumber(GeneratedMeshUserNumber, NumberOfElements, &
    & CMISSGeneratedMeshNumberOfElementsSetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshNumberOfElementsSetNumberC

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshNumberOfElementsSetPtrC(GeneratedMeshPtr,NumberOfElements) BIND(C, NAME = &
  & "CMISSGeneratedMeshNumberOfElementsSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to set the number of elements to.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements(:) !<NumberOfElements(i). On return, the number of elements in the i'th dimension of the generated mesh to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshNumberOfElementsSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshNumberOfElementsSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshNumberOfElementsSetObj(GeneratedMesh, NumberOfElements, CMISSGeneratedMeshNumberOfElementsSetPtrC)
      ELSE
        CMISSGeneratedMeshNumberOfElementsSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshNumberOfElementsSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshNumberOfElementsSetPtrC

  !
  !================================================================================================================================
  !

  !>Returns the origin of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshOriginGetNumberC(GeneratedMeshUserNumber,Origin) BIND(C, NAME = "CMISSGeneratedMeshOriginGetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to get the origin for, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Origin(:) !<Origin(i). On return, the origin of the i'th dimension of the generated mesh to get for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshOriginGetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshOriginGetNumber(GeneratedMeshUserNumber, Origin, CMISSGeneratedMeshOriginGetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshOriginGetNumberC

  !
  !================================================================================================================================
  !

  !>Returns the origin of a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshOriginGetPtrC(GeneratedMeshPtr,Origin) BIND(C, NAME = "CMISSGeneratedMeshOriginGetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to set the number of elements to.
    REAL(C_DOUBLE), INTENT(OUT) :: Origin(:) !<Origin(i). On return, the origin of the i'th dimension of the generated mesh to get for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshOriginGetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshOriginGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshOriginGetObj(GeneratedMesh, Origin, CMISSGeneratedMeshOriginGetPtrC)
      ELSE
        CMISSGeneratedMeshOriginGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshOriginGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshOriginGetPtrC

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshOriginSetNumberC(GeneratedMeshUserNumber,Origin) BIND(C, NAME = "CMISSGeneratedMeshOriginSetNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the origin to, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Origin(:) !<Origin(i). On return, the origin of the i'th dimension of the generated mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshOriginSetNumberC !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshOriginSetNumber(GeneratedMeshUserNumber, Origin, CMISSGeneratedMeshOriginSetNumberC)

    RETURN

  END FUNCTION CMISSGeneratedMeshOriginSetNumberC

    !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshOriginSetPtrC(GeneratedMeshPtr,Origin) BIND(C, NAME = "CMISSGeneratedMeshOriginSetObj")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to set the number of elements to.
    REAL(C_DOUBLE), INTENT(OUT) :: Origin(:) !<Origin(i). On return, the origin of the i'th dimension of the generated mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshOriginSetPtrC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshOriginSetPtrC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshOriginSetObj(GeneratedMesh, Origin, CMISSGeneratedMeshOriginSetPtrC)
      ELSE
        CMISSGeneratedMeshOriginSetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshOriginSetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshOriginSetPtrC

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
    TYPE(CMISSMeshElementsType), POINTER :: Mesh
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
