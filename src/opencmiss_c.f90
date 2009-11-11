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

  PUBLIC CMISSCoordinateSystemTypeFinaliseC,CMISSCoordinateSystemTypeInitialiseC

  PUBLIC CMISSRegionTypeFinaliseC,CMISSRegionTypeInitialiseC
  
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

  !>Finalises a CMISSCoordinateSystemType object for C.
  FUNCTION CMISSCoordinateSystemTypeFinaliseC(CoordinateSystemTypePtr)  BIND(C,NAME="CMISSCoordinateSystemTypeFinalise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: CoordinateSystemTypePtr
    !Function variable
    INTEGER(C_INT) :: CMISSCoordinateSystemTypeFinaliseC
    !Local variables
    TYPE(CMISSCoordinateSystemType), POINTER :: CoordinateSystemType
    
    CMISSCoordinateSystemTypeFinaliseC=CMISSNoError
    IF(C_ASSOCIATED(CoordinateSystemTypePtr)) THEN
      CALL C_F_POINTER(CoordinateSystemTypePtr,CoordinateSystemType)
      IF(ASSOCIATED(CoordinateSystemType)) THEN
        CALL CMISSCoordinateSystemTypeFinalise(CoordinateSystemType,CMISSCoordinateSystemTypeFinaliseC)
        DEALLOCATE(CoordinateSystemType)
        CoordinateSystemTypePtr=C_NULL_PTR
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSCoordinateSystemTypeFinaliseC
 
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSCoordinateSystemType object for C.
  FUNCTION CMISSCoordinateSystemTypeInitialiseC(CoordinateSystemTypePtr)  BIND(C,NAME="CMISSCoordinateSystemTypeInitialise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: CoordinateSystemTypePtr
    !Function variable
    INTEGER(C_INT) :: CMISSCoordinateSystemTypeInitialiseC
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

  !>Finalises a CMISSRegionType object for C.
  FUNCTION CMISSRegionTypeFinaliseC(RegionTypePtr) BIND(C,NAME="CMISSRegionTypeFinalise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: RegionTypePtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionTypeFinaliseC
    !Local variables
    TYPE(CMISSRegionType), POINTER :: RegionType
    
    CMISSRegionTypeFinaliseC=CMISSNoError
    IF(C_ASSOCIATED(RegionTypePtr)) THEN
      CALL C_F_POINTER(RegionTypePtr,RegionType)
      IF(ASSOCIATED(RegionType)) THEN
        CALL CMISSRegionTypeFinalise(RegionType,CMISSRegionTypeFinaliseC)
        DEALLOCATE(RegionType)
        RegionTypePtr=C_NULL_PTR
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionTypeFinaliseC
 
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSRegionType object for C.
  FUNCTION CMISSRegionTypeInitialiseC(RegionTypePtr) BIND(C,NAME="CMISSRegionTypeInitialise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: RegionTypePtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionTypeInitialiseC
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
!! REGION_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the process of creating a region identified by a pointer for C.
  FUNCTION CMISSRegionCreateFinishCPtr(RegionPtr) BIND(C,NAME="CMISSRegionCreateFinish")
  
    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateFinishCPtr
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region

    CMISSRegionCreateFinishCPtr=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSRegionCreateFinish(Region,CMISSRegionCreateFinishCPtr)
      ELSE
        CMISSRegionCreateFinishCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionCreateFinishCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionCreateFinishCPtr

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
  FUNCTION CMISSRegionCreateStartCPtr(RegionUserNumber,ParentRegionPtr,RegionPtr) BIND(C,NAME="CMISSRegionCreateStart")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    TYPE(C_PTR), VALUE, INTENT(IN) :: ParentRegionPtr
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateStartCPtr
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region,ParentRegion

    CMISSRegionCreateStartCPtr=CMISSNoError
    IF(C_ASSOCIATED(ParentRegionPtr)) THEN
      CALL C_F_POINTER(ParentRegionPtr,ParentRegion)
      IF(ASSOCIATED(ParentRegion)) THEN        
        IF(C_ASSOCIATED(RegionPtr)) THEN
          CALL C_F_POINTER(RegionPtr,Region)
          IF(ASSOCIATED(Region)) THEN        
            CALL CMISSRegionCreateStart(RegionUserNumber,ParentRegion,Region,CMISSRegionCreateStartCPtr)
          ELSE
            CMISSRegionCreateStartCPtr=CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSRegionCreateStartCPtr=CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSRegionCreateStartCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionCreateStartCPtr=CMISSPointerIsNULL
    ENDIF
    
    RETURN
    
  END FUNCTION CMISSRegionCreateStartCPtr

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
