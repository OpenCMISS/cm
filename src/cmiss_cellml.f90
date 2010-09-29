!> \file
!> $Id$
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief This module is a OpenCMISS(cm) buffer module to OpenCMISS(cellml).
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

!> This module is a OpenCMISS(cm) buffer module to OpenCMISS(cellml).
MODULE CMISS_CELLML

  USE ISO_C_BINDING

#ifdef USECELLML
  USE CELLML_MODEL_DEFINITION
#endif

  !Module imports
  USE BASE_ROUTINES
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE INPUT_OUTPUT
  USE KINDS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup CELLML_FieldTypes CMISS_CELLML::FieldTypes
  !> \brief CellML field parameter types
  !> \see CMISS_CELLML,OPENCMISS_CellMLFieldTypes
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD_TYPE = 1 !<CellML models field type \see CELLML_FieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_STATE_FIELD_TYPE = 2 !<CellML state field type \see CELLML_FieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_INTERMEDIATE_FIELD_TYPE = 3 !<CellML intermediate field type \see CELLML_FieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_PARAMETERS_FIELD_TYPE = 4 !<CellML parameters field type \see CELLML_FieldTypes,CMISS_CELLML
  !>@}

  !Module types

  !Module variables

  TYPE(CELLML_ENVIRONMENTS_TYPE), TARGET :: CELLML_ENVIRONMENTS
 
  !Interfaces

  INTERFACE CELLML_MODEL_IMPORT
    MODULE PROCEDURE CELLML_MODEL_IMPORT_C
    MODULE PROCEDURE CELLML_MODEL_IMPORT_VS
  END INTERFACE !CELLML_MODEL_IMPORT
  
  INTERFACE CELLML_CREATE_FIELD_TO_CELLML_MAP
    MODULE PROCEDURE CELLML_CREATE_FIELD_TO_CELLML_MAP_C
    MODULE PROCEDURE CELLML_CREATE_FIELD_TO_CELLML_MAP_VS
  END INTERFACE !CELLML_CREATE_FIELD_TO_CELLML_MAP

  INTERFACE CELLML_FIELD_COMPONENT_GET
    MODULE PROCEDURE CELLML_FIELD_COMPONENT_GET_C
    MODULE PROCEDURE CELLML_FIELD_COMPONENT_GET_VS
  END INTERFACE !CELLML_FIELD_COMPONENT_GET

  INTERFACE CELLML_INTERMEDIATE_FIELD_ADD
    MODULE PROCEDURE CELLML_INTERMEDIATE_FIELD_ADD_C
    MODULE PROCEDURE CELLML_INTERMEDIATE_FIELD_ADD_VS
  END INTERFACE !CELLML_INTERMEDIATE_FIELD_ADD

  INTERFACE CELLML_PARAMETER_ADD
    MODULE PROCEDURE CELLML_PARAMETER_ADD_C
    MODULE PROCEDURE CELLML_PARAMETER_ADD_VS
  END INTERFACE !CELLML_PARAMETER_ADD

  PUBLIC CELLML_MODELS_FIELD_TYPE,CELLML_STATE_FIELD_TYPE,CELLML_INTERMEDIATE_FIELD_TYPE,CELLML_PARAMETERS_FIELD_TYPE
  
  PUBLIC CELLML_CREATE_START,CELLML_CREATE_FINISH

  PUBLIC CELLML_DESTROY
  
  PUBLIC CELLML_MODEL_IMPORT

  PUBLIC CELLML_CREATE_FIELD_TO_CELLML_MAP

  PUBLIC CELLML_MODELS_FIELD_CREATE_START,CELLML_MODELS_FIELD_CREATE_FINISH,CELLML_MODELS_FIELD_GET

  PUBLIC CELLML_STATE_FIELD_CREATE_START,CELLML_STATE_FIELD_CREATE_FINISH,CELLML_STATE_FIELD_GET

  PUBLIC CELLML_FIELD_COMPONENT_GET

  PUBLIC CELLML_INTERMEDIATE_FIELD_ADD,CELLML_INTERMEDIATE_FIELD_CREATE_FINISH,CELLML_INTERMEDIATE_FIELD_CREATE_START, &
    & CELLML_INTERMEDIATE_FIELD_GET

  PUBLIC CELLML_PARAMETER_ADD,CELLML_PARAMETERS_CREATE_FINISH,CELLML_PARAMETERS_CREATE_START

  PUBLIC CELLML_PARAMETERS_FIELD_CREATE_START,CELLML_PARAMETERS_FIELD_CREATE_FINISH,CELLML_PARAMETERS_FIELD_GET

  PUBLIC CELLML_GENERATE

  PUBLIC CELLML_USER_NUMBER_FIND,CELLML_MODEL_USER_NUMBER_FIND

  PUBLIC CELLML_ENVIRONMENTS_FINALISE,CELLML_ENVIRONMENTS_INITIALISE
  
CONTAINS

  !
  !=================================================================================================================================
  !

  !> Destroys the given CellML environment.
  SUBROUTINE CELLML_DESTROY(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_DESTROY",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CALL CELLML_FINALISE(CELLML,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("CellML is not associated.",ERR,ERROR,*999)
    END IF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_DESTROY")
    RETURN
999 CALL ERRORS("CELLML_DESTROY",ERR,ERROR)
    CALL EXITS("CELLML_DESTROY")
    RETURN 1
  END SUBROUTINE CELLML_DESTROY

  !
  !=================================================================================================================================
  !

  !>Set up the CellML environment in the given field.
  !! For a given region, create a CellML environment that will be used to define CellML models in openCMISS. This will simply create and initialise an empty CellML environment object in the specified region with the specified unique identifier number. Also set some flag to indicate the CellML environment object is in the process of being created and should not yet be used.
  !! \todo Should be passing in a region? or leave that for the field mapping?
  SUBROUTINE CELLML_CREATE_START(CELLML_USER_NUMBER,REGION,CELLML,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CELLML_USER_NUMBER !<The user specified ID for this CellML environment.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: REGION !<The region this CellML environment belongs to.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The newly created CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: cellml_idx
    TYPE(CELLML_TYPE), POINTER :: NEW_CELLML
    TYPE(CELLML_PTR_TYPE), POINTER :: NEW_CELLML_ENVIRONMENTS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_CELLML)

    CALL ENTERS("CELLML_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

  IF(.NOT. ASSOCIATED(REGION)) THEN
    CALL FLAG_ERROR("Region for a new CellML environment must be defined.",ERR,ERROR,*999)
  ENDIF
  IF(ASSOCIATED(CELLML)) THEN
    CALL FLAG_ERROR("CellML is already associated.",ERR,ERROR,*999)
  ELSE
    NULLIFY(CELLML)
    CALL CELLML_USER_NUMBER_FIND(CELLML_USER_NUMBER,CELLML,ERR,ERROR,*999)
    IF(ASSOCIATED(CELLML)) THEN
      LOCAL_ERROR="CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML_USER_NUMBER,"*",ERR,ERROR))// &
        & " has already been created."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ELSE
      !Allocate new cellml
      ALLOCATE(NEW_CELLML,STAT=ERR)
      IF (ERR/=0) CALL FLAG_ERROR("Could not allocate new CellML environment.",ERR,ERROR,*999)
      !Initialise cellml
      CALL CELLML_INITIALISE(NEW_CELLML,ERR,ERROR,*999)
      !Set cellml defaults
      NEW_CELLML%USER_NUMBER = CELLML_USER_NUMBER
      NEW_CELLML%GLOBAL_NUMBER = CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1
      NEW_CELLML%ENVIRONMENTS => CELLML_ENVIRONMENTS
      NEW_CELLML%REGION => REGION
      NEW_CELLML%CELLML_FINISHED = .FALSE.
      NULLIFY(NEW_CELLML%MODELS)
      ALLOCATE(NEW_CELLML%MODELS,STAT=ERR)
      IF (ERR/=0) CALL FLAG_ERROR("CellML environment models could not be allocated.",ERR,ERROR,*999)
      NEW_CELLML%MODELS%NUMBER_OF_MODELS = 0
      NULLIFY(NEW_CELLML%MODELS%MODELS)
      NEW_CELLML%MODELS%CELLML => NEW_CELLML
      !Add cellml to the list of cellml environments
      ALLOCATE(NEW_CELLML_ENVIRONMENTS(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new CellML environments.",ERR,ERROR,*999)
      DO cellml_idx=1,CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS
        NEW_CELLML_ENVIRONMENTS(cellml_idx)%PTR=>CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR
      ENDDO !cellml_idx
      NEW_CELLML_ENVIRONMENTS(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1)%PTR=>NEW_CELLML
      IF(ASSOCIATED(CELLML_ENVIRONMENTS%ENVIRONMENTS)) DEALLOCATE(CELLML_ENVIRONMENTS%ENVIRONMENTS)
      CELLML_ENVIRONMENTS%ENVIRONMENTS=>NEW_CELLML_ENVIRONMENTS
      CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS=CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1
      CELLML=>NEW_CELLML
    ENDIF
  ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_CELLML)) CALL CELLML_DESTROY(NEW_CELLML,ERR,ERROR,*998)
998 NULLIFY(CELLML)
    CALL ERRORS("CELLML_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_START

  !
  !=================================================================================================================================
  !

  !> Finish creating the CellML environment.
  !! Check the provided CellML environment object and if it all looks good clear the "in progress" flag to indicate the object is now ready for use.
  SUBROUTINE CELLML_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object to check and finialise creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL ENTERS("CELLML_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
!      IF(ASSOCIATED(CELLML%SOURCE_FIELD)) THEN
!        !Insert this CellML environment into the source field, deleting any existing CellML environment
!        IF(ASSOCIATED(CELLML%SOURCE_FIELD%CELLML)) CALL CELLML_DESTROY(CELLML%SOURCE_FIELD%CELLML,ERR,ERROR,*999)
!        CELLML%SOURCE_FIELD%CELLML => CELLML
        CELLML%CELLML_FINISHED = .TRUE.
!      ELSE
!        CALL FLAG_ERROR("CellML source field not associated.",ERR,ERROR,*999)
!      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML environment and deallocate all memory.
  SUBROUTINE CELLML_FINALISE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
!      IF(ASSOCIATED(CELLML%SOURCE_FIELD)) THEN
!        IF(ASSOCIATED(CELLML%SOURCE_FIELD%CELLML)) NULLIFY(CELLML%SOURCE_FIELD%CELLML)
!      ENDIF
      IF(ASSOCIATED(CELLML%MODELS)) CALL CELLML_MODELS_FINALISE(CELLML%MODELS,ERR,ERROR,*999)
      DEALLOCATE(CELLML)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_FINALISE

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML environment and deallocate all memory.
  SUBROUTINE CELLML_INITIALISE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    CALL ENTERS("CELLML_INITIALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
        NULLIFY(CELLML%REGION)
        NULLIFY(CELLML%MODELS_FIELD)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INITIALISE")
    RETURN
999 CALL ERRORS("CELLML_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML models deallocate all memory.
  SUBROUTINE CELLML_MODELS_FINALISE(CELLML_MODELS,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_MODELS_TYPE), POINTER :: CELLML_MODELS !<A pointer to the CellML models to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_MODELS_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODELS)) THEN

    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODELS_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FINALISE

  !
  !=================================================================================================================================
  !

  !>Import the specified CellML model into the given CellML environment object.
  !! Here we load specified CellML models into the CellML environment object. Will be called for each model required for use with this CellML environment.
  !! - should do some level of validation when the model is loaded
  !! - see URI notes below...
  SUBROUTINE CELLML_MODEL_IMPORT_C(MODEL_USER_NUMBER,CELLML,URI,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The unique identifier for this model within the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object into which we want to import the specified model.
    CHARACTER(LEN=*) :: URI !<The (absolute? relative?) URI of the model to import. As per tracker item 2013 comment 8 the URI should now simply point to a CellML document. Can use a relative URL which will be interpreted relative to the CWD of the executed application.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: cellml_idx
    TYPE(CELLML_MODEL_TYPE), POINTER :: NEW_MODEL
    TYPE(CELLML_MODEL_PTR_TYPE), POINTER :: NEW_MODELS(:)
    !INTEGER (C_INT) CODE
    !TYPE(C_PTR) :: CELLML_MODEL
    CHARACTER(256) :: C_URI
    INTEGER(INTG) :: C_URI_L
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_MODEL_IMPORT_C",ERR,ERROR,*999)

#ifdef USECELLML

    NULLIFY(NEW_MODEL)

    IF (ASSOCIATED(CELLML)) THEN
        CALL CELLML_MODEL_USER_NUMBER_FIND(MODEL_USER_NUMBER,CELLML,NEW_MODEL,ERR,ERROR,*999)
        IF(ASSOCIATED(NEW_MODEL)) THEN
            LOCAL_ERROR="CellML model number "//TRIM(NUMBER_TO_VSTRING(MODEL_USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ELSE
            !Allocate new CellML model
            ALLOCATE(NEW_MODEL,STAT=ERR)
            IF (ERR/=0) CALL FLAG_ERROR("Could not allocate new CellML model.",ERR,ERROR,*999)
            !set up new model object
            C_URI_L = LEN_TRIM(URI)
            WRITE(C_URI,'(A,A)') URI(1:C_URI_L),C_NULL_CHAR
            NEW_MODEL%PTR = CREATE_CELLML_MODEL_DEFINITION(C_URI)
            NEW_MODEL%USER_NUMBER = MODEL_USER_NUMBER
            NEW_MODEL%GLOBAL_NUMBER = CELLML%MODELS%NUMBER_OF_MODELS + 1
            ! Add model to this CellML environment's list of models
            ALLOCATE(NEW_MODELS(CELLML%MODELS%NUMBER_OF_MODELS+1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new CellML models array.",ERR,ERROR,*999)
            DO cellml_idx=1,CELLML%MODELS%NUMBER_OF_MODELS
                NEW_MODELS(cellml_idx)%PTR=>CELLML%MODELS%MODELS(cellml_idx)%PTR
            ENDDO !cellml_idx
            NEW_MODELS(CELLML%MODELS%NUMBER_OF_MODELS+1)%PTR=>NEW_MODEL
            IF(ASSOCIATED(CELLML%MODELS%MODELS)) DEALLOCATE(CELLML%MODELS%MODELS)
            CELLML%MODELS%MODELS => NEW_MODELS
            CELLML%MODELS%NUMBER_OF_MODELS = CELLML%MODELS%NUMBER_OF_MODELS+1
        ENDIF
    ELSE
        CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODEL_IMPORT_C")
    RETURN
999 CALL ERRORS("CELLML_MODEL_IMPORT_C",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_IMPORT_C")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_IMPORT_C

  !
  !=================================================================================================================================
  !

  !>Import the specified CellML model into the given CellML environment object.
  !! Here we load specified CellML models into the CellML environment object. Will be called for each model required for use with this CellML environment.
  !! - should do some level of validation when the model is loaded
  !! - see URI notes below...
  SUBROUTINE CELLML_MODEL_IMPORT_VS(MODEL_USER_NUMBER,CELLML,URI,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The unique identifier for this model within the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object into which we want to import the specified model.
    TYPE(VARYING_STRING), INTENT(IN) :: URI !<The (absolute? relative?) URI of the model to import. As per tracker item 2013 comment 8 the URI should now simply point to a CellML document. Can use a relative URL which will be interpreted relative to the CWD of the executed application.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODEL_IMPORT_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_MODEL_IMPORT(MODEL_USER_NUMBER,CELLML,CHAR(URI),ERR,ERROR,*999)

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODEL_IMPORT_VS")
    RETURN
999 CALL ERRORS("CELLML_MODEL_IMPORT_VS",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_IMPORT_VS")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_IMPORT_VS

  !
  !=================================================================================================================================
  !

  !>Create a field variable component to CellML model variable map.
  SUBROUTINE CELLML_CREATE_FIELD_TO_CELLML_MAP_C(CELLML,FIELD,FIELD_VARIABLE,FIELD_COMPONENT,MODEL_USER_NUMBER,VARIABLE_ID, &
  & ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field to map from.
    INTEGER(INTG), INTENT(IN) :: FIELD_VARIABLE !<The field variable to map from.
    INTEGER(INTG), INTENT(IN) :: FIELD_COMPONENT !<The component number to map from the given field variable.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the CellML model to map to.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The ID of the CellML variable in the given model to map to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_CREATE_FIELD_TO_CELLML_MAP_C",ERR,ERROR,*999)

#ifdef USECELLML

    write(*,*) 'Mapping variable: ',FIELD_VARIABLE,' to model number: ',MODEL_USER_NUMBER,' with ID: ',VARIABLE_ID

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_FIELD_TO_CELLML_MAP_C")
    RETURN
999 CALL ERRORS("CELLML_CREATE_FIELD_TO_CELLML_MAP_C",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_FIELD_TO_CELLML_MAP_C")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_FIELD_TO_CELLML_MAP_C

  !
  !=================================================================================================================================
  !

  !>Create a field variable component to CellML model variable map.
  SUBROUTINE CELLML_CREATE_FIELD_TO_CELLML_MAP_VS(CELLML,FIELD,FIELD_VARIABLE,FIELD_COMPONENT,MODEL_USER_NUMBER,VARIABLE_ID, &
  & ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field to map from.
    INTEGER(INTG), INTENT(IN) :: FIELD_VARIABLE !<The field variable to map from.
    INTEGER(INTG), INTENT(IN) :: FIELD_COMPONENT !<The component number to map from the given field variable.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the CellML model to map to.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_ID !<The ID of the CellML variable in the given model to map to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_CREATE_FIELD_TO_CELLML_MAP_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_CREATE_FIELD_TO_CELLML_MAP(CELLML,FIELD,FIELD_VARIABLE,FIELD_COMPONENT,MODEL_USER_NUMBER,CHAR(VARIABLE_ID), &
    & ERR,ERROR,*999)

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_FIELD_TO_CELLML_MAP_VS")
    RETURN
999 CALL ERRORS("CELLML_CREATE_FIELD_TO_CELLML_MAP_VS",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_FIELD_TO_CELLML_MAP_VS")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_FIELD_TO_CELLML_MAP_VS

  !
  !=================================================================================================================================
  !

  !>Start the creation of the models field for the given CellML environment.
  !! This will create the models field for the given CellML environment object. The models field is used to associate specific models defined for this CellML environment with each of the degrees of freedom for the (source) field for which this CellML environment is defined.
  !! - what to do if models field already exists? exists but has a different user number?
  !! - perhaps a better way to do this is through the mapping of CellML model variables to field components. This might provide more flexibility, including the use of multiple models for different field components at a single DOF. Otherwise, still need some way to map model variables to field components.
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_START(MODEL_FIELD_USER_NUMBER,CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: MODEL_FIELD_USER_NUMBER !<The unique identifier for the models field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to create the models field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On Return, the created models field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR,MESH_COMPONENT
    TYPE(REGION_TYPE), POINTER :: REGION,MODELS_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML_NEEDS_REDOING

    IF(ASSOCIATED(CELLML) .AND. ASSOCIATED(CELLML%SOURCE_FIELD)) THEN
      IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
        CALL FLAG_ERROR("The CellML environment models field is already associated",ERR,ERROR,*999)
      ELSE
        REGION=>CELLML%SOURCE_FIELD%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(FIELD)) THEN
            !Check the materials field has been finished
            IF(FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(MODEL_FIELD_USER_NUMBER/=FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified models field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(MODEL_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified models field of "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              MODELS_FIELD_REGION=>FIELD%REGION
              IF(ASSOCIATED(MODELS_FIELD_REGION)) THEN
                !Check the field is defined on the same region as the source field
                IF(MODELS_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified models field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(MODELS_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified source field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified models field has the same decomposition as the source field
                  IF(.NOT.ASSOCIATED(CELLML%SOURCE_FIELD%DECOMPOSITION,FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified models field does not have the same decomposition as the source "// &
                      & "field for the specified CellML environment.",ERR,ERROR,*999)
                  ENDIF
              ELSE
                CALL FLAG_ERROR("The specified models field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified models field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(MODEL_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified models field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(MODEL_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Initialise the environment models field
          IF(.NOT. ASSOCIATED(FIELD)) THEN
            ! need to create a new models field from scratch
            NULLIFY(FIELD)
            CALL FIELD_CREATE_START(MODEL_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            CALL FIELD_DATA_TYPE_SET_AND_LOCK(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,ERR,ERROR,*999)
            CALL FIELD_LABEL_SET_AND_LOCK(FIELD,"CellMLModelsField",ERR,ERROR,*999)
            CALL FIELD_TYPE_SET_AND_LOCK(FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
            CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(FIELD,CELLML%SOURCE_FIELD%DECOMPOSITION,ERR,ERROR,*999)
            CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(FIELD,CELLML%SOURCE_FIELD%GEOMETRIC_FIELD,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(FIELD,1,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_LABEL_SET_AND_LOCK(FIELD,FIELD_U_VARIABLE_TYPE,"ModelMap",ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
            CALL FIELD_COMPONENT_LABEL_SET_AND_LOCK(FIELD,FIELD_U_VARIABLE_TYPE,1,"ModelUserNumber",ERR,ERROR,*999)
            CALL FIELD_COMPONENT_MESH_COMPONENT_GET(CELLML%SOURCE_FIELD, FIELD_U_VARIABLE_TYPE, 1, MESH_COMPONENT, ERR, ERROR, *999)
            CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(FIELD,FIELD_U_VARIABLE_TYPE,1,1,ERR,ERROR,*999)
            CALL FIELD_CREATE_FINISH(FIELD,ERR,ERROR,*999)
          ENDIF
          ! save the specified models field
          CELLML%MODELS_FIELD => FIELD
        ELSE
          CALL FLAG_ERROR("CellML environment source field region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated",ERR,ERROR,*999)
    ENDIF

#else
    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
#endif
    
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finish the creation of the models field for the given CellML environment.
  !! This will finalise the creation of the models field for the given CellML environment object.
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to finish creation of the models field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

    ! nothing to do?

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODELS_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Fetch the models field for the given CellML environment.
  !! Check the models field is correctly defined and return it for the user to define. The user will be able to use the returned field object to assign models to DOFs.
  SUBROUTINE CELLML_MODELS_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the models field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the models field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODELS_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_GET

  !
  !=================================================================================================================================
  !

  !>Start the creation of the state field for the given CellML environment.
  !! We want to create a single field which will be used to store the state variables of each of the models in the given CellML environment at each of the DOFs for the environment's base field.
  !! - one field for all models with a component in the field for each variable (perhaps simply use the largest model to set the number of components for all DOFs and some will not be used).
  !! - can a single field have different components at different DOFs (grid points)?
  !! - perhaps need a state field for each model and each field is only defined at the DOFs valid for that model?
  !! - this has to be called after the models field is defined.
  SUBROUTINE CELLML_STATE_FIELD_CREATE_START(STATE_FIELD_USER_NUMBER,CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: STATE_FIELD_USER_NUMBER !<The unique identifier for the state field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to create the state field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On Return, the created state field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_STATE_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finialse the creation of the state field for the given CellML environment.
  !! Finish creating the state variable field for the provided CellML environment.
  !! - default field values are set for all components in the field at all DOFs based on the initial_value's from the CellML model
  SUBROUTINE CELLML_STATE_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the creation of the state field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_STATE_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Fetch the state field for the given CellML environment.
  !! Check the state field is correctly defined and return it for the user to access - either for getting the current values or to set initial conditions.
  !! - need some way to get from variable URIs in the CellML models to components in the returned state field.
  !! - standard field routines should be used to get/set field component values.
  SUBROUTINE CELLML_STATE_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the state field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the state field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_STATE_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_GET

  !
  !=================================================================================================================================
  !

  !>Find the component ID in the given field for the variable defined by the given URI in the provided CellML environment.
  !! This generic routine will be used to map variable URI's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model URI also needed?
  SUBROUTINE CELLML_FIELD_COMPONENT_GET_C(CELLML,CELLML_FIELD_TYPE,URI,COMPONENT_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CMISS_CELLML
    CHARACTER(LEN=*), INTENT(IN) :: URI !<The URI of the model variable which needs to be located in the provided field.
    INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On return, the field component for the model variable defined by the given URI.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_FIELD_COMPONENT_GET_C",ERR,ERROR,*999)

#ifdef USECELLML

    COMPONENT_USER_NUMBER=0

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FIELD_COMPONENT_GET_C")
    RETURN
999 CALL ERRORS("CELLML_FIELD_COMPONENT_GET_C",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_COMPONENT_GET_C")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_COMPONENT_GET_C

  !
  !=================================================================================================================================
  !

  !>Find the component ID in the given field for the variable defined by the given URI in the provided CellML environment.
  !! This generic routine will be used to map variable URI's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model URI also needed?
  SUBROUTINE CELLML_FIELD_COMPONENT_GET_VS(CELLML,CELLML_FIELD_TYPE,URI,COMPONENT_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CMISS_CELLML
    TYPE(VARYING_STRING), INTENT(IN) :: URI !<The URI of the model variable which needs to be located in the provided field.
    INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On return, the field component for the model variable defined by the given URI.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_FIELD_COMPONENT_GET_VS",ERR,ERROR,*999)

#ifdef USECELLML

    COMPONENT_USER_NUMBER=0

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FIELD_COMPONENT_GET_VS")
    RETURN
999 CALL ERRORS("CELLML_FIELD_COMPONENT_GET_VS",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_COMPONENT_GET_VS")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_COMPONENT_GET_VS

  !
  !=================================================================================================================================
  !

  !>Create a field used to store intermediate variables of interest.
  !! The intermediate field created for this environment provides access to the value of nominated intermediate variables from the models available in this CellML environment.
  !! - this field will be similar in structure to the state field, but only contains the nominated variables.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_START(INTERMEDIATE_FIELD_USER_NUMBER,CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: INTERMEDIATE_FIELD_USER_NUMBER !<The unique identifier for the intermediate field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On return, the created intermediates field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finialse the creation of the intermediate field for the given CellML environment.
  !! Finish creating the intermediate variable field for the provided CellML environment.
  !! - check for valid variables being defined?
  !! - maybe delete the field if no valid variables/components exist?
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the creation of the intermediate field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Add a specific variable to this CellML environment's intermediate field.
  !! Nominate a specific variable from a model in this environment for inclusion in the environment's intermediate field. This will ensure the value of this variable is available at each of the environment's (source) field DOF's for which the specified model is valid.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD_C(CELLML,MODEL_USER_NUMBER,VARIABLE_URI,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to add the nominated variable to the intermediate field.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the model (simulation?) in which to look for the variable being nominated.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_URI !<The URI of the variable in the specified model to include in the intermediate field for this CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_ADD_C",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD_C")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_ADD_C",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD_C")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD_C

  !
  !=================================================================================================================================
  !

  !>Add a specific variable to this CellML environment's intermediate field.
  !! Nominate a specific variable from a model in this environment for inclusion in the environment's intermediate field. This will ensure the value of this variable is available at each of the environment's (source) field DOF's for which the specified model is valid.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD_VS(CELLML,MODEL_USER_NUMBER,VARIABLE_URI,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to add the nominated variable to the intermediate field.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the model (simulation?) in which to look for the variable being nominated.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_URI !<The URI of the variable in the specified model to include in the intermediate field for this CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_ADD_VS",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD_VS")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_ADD_VS",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD_VS")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD_VS

  !
  !=================================================================================================================================
  !

  !>Fetch the intermediate field for the given CellML environment.
  !! Check the intermediate field is correctly defined and return it for the user to access.
  !! - is there a way to ensure the field returned to the calling routine is "read-only"?
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the intermediate field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the intermediate field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_GET

  !
  !=================================================================================================================================
  !

  !>Start the parameters definition process.
  !! Initialise the CellML environment ready for the nomination of parameters from this environment's models that will be overridden by field values.
  !! - what to do if the paramters field has already be defined?
  SUBROUTINE CELLML_PARAMETERS_CREATE_START(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we will be defining parameter overrides.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETERS_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finialse the parameters definition process.
  !! Indicates that the user has added all parameter overrides and allows the CellML environment to be further processed.
  SUBROUTINE CELLML_PARAMETERS_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the parameter overrides.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETERS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Nominate a specific parameter in CellML environment to override.
  !! Nominate a specific parameter variable from a model in this environment which will have its value overridden by the envronment's parameter field.
  SUBROUTINE CELLML_PARAMETER_ADD_C(CELLML,MODEL_USER_NUMBER,VARIABLE_URI,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the model (simulation?) in which to look for the variable being nominated.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_URI !<The URI of the variable in the specified model to specify for overriding.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETER_ADD_C",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETER_ADD_C")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_ADD_C",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_ADD_C")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETER_ADD_C

  !
  !=================================================================================================================================
  !

  !>Nominate a specific parameter in CellML environment to override.
  !! Nominate a specific parameter variable from a model in this environment which will have its value overridden by the envronment's parameter field.
  SUBROUTINE CELLML_PARAMETER_ADD_VS(CELLML,MODEL_USER_NUMBER,VARIABLE_URI,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the model (simulation?) in which to look for the variable being nominated.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_URI !<The URI of the variable in the specified model to specify for overriding.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETER_ADD_VS",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETER_ADD_VS")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_ADD_VS",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_ADD_VS")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETER_ADD_VS

  !
  !=================================================================================================================================
  !

  !>Create the parameters field for this CellML environment.
  !! Check that parameters have been nominated and create the environment's user accessible parameters field. The user accessible parameters field is used by the user to define parameter overrides and is distinct from the internally used parameters field which is always defined for the environment's (source) field DOF's.
  !! - similar to the models field.
  !! - should initialise all components of the field as constant fields with a value as specified in the corresponding CellML models imported into the environment?
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_START(PARAMETERS_FIELD_USER_NUMBER,CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARAMETERS_FIELD_USER_NUMBER !<The unique identifier for the parameters field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we will be defining the parameters field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On return, the created parameters field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finialse the parameters field definition.
  !! - here is where we sync the user level parameters field to the internal parameters representation?
  !! - or maybe we simply create the internal field?
  !! - the internal field will contain components for all parameters in all models?
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the parameters field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Fetch the user parameters field for the given CellML environment.
  !! Check the parameters field is correctly defined and return it for the user to access in order to define the desired parameter overrides.
  SUBROUTINE CELLML_PARAMETERS_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the prameters field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the parameters field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETERS_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_GET

  !
  !=================================================================================================================================
  !

  !>Validate and instantiate the specified CellML environment.
  !! Users should call this routine once they have set up the CellML environment to allow the CellML environment to be validated and instantiated into computable code.
  !! - generate and compile code.
  !! - distribute things amonst processors?
  !! - create the internal fields?
  !! - does this method need to get called or can it be done implicitly as needed by other solution procedures?
  SUBROUTINE CELLML_GENERATE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to generate computable code.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_GENERATE",ERR,ERROR,*999)

#ifdef USECELLML

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_GENERATE")
    RETURN
999 CALL ERRORS("CELLML_GENERATE",ERR,ERROR)
    CALL EXITS("CELLML_GENERATE")
    RETURN 1
  END SUBROUTINE CELLML_GENERATE

  !
  !=================================================================================================================================
  !

  !>Finds and returns in CELLML a pointer to the CellML environment identified by USER_NUMBER. If no CellML environment with that USER_NUMBER exists CELLML is left nullified.
  SUBROUTINE CELLML_USER_NUMBER_FIND(USER_NUMBER,CELLML,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<On return a pointer to the CellML environment with the given user number. If no CellML environment with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: cellml_idx

    CALL ENTERS("CELLML_USER_NUMBER_FIND",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CALL FLAG_ERROR("CellML is already associated.",ERR,ERROR,*999)
    ELSE
      cellml_idx=1
      DO WHILE(cellml_idx<=CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS.AND..NOT.ASSOCIATED(CELLML))
        IF(CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          CELLML=>CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR
        ELSE
          cellml_idx=cellml_idx+1
        ENDIF
      ENDDO
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("CELLML_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("CELLML_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE CELLML_USER_NUMBER_FIND

  !>Finds and returns in CELLML a pointer to the CellML model identified by USER_NUMBER. If no CellML environment with that USER_NUMBER exists CELLML is left nullified.
  SUBROUTINE CELLML_MODEL_USER_NUMBER_FIND(USER_NUMBER,CELLML,CELLML_MODEL,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment in which to search for the specified user number.
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL !<On return a pointer to the CellML model with the given user number within the given CellML environment. If no CellML model with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: cellml_idx

    CALL ENTERS("CELLML_MODEL_USER_NUMBER_FIND",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL)) THEN
      CALL FLAG_ERROR("CellML model is already associated.",ERR,ERROR,*999)
    ELSE
      cellml_idx=1
      DO WHILE(cellml_idx<=CELLML%MODELS%NUMBER_OF_MODELS.AND..NOT.ASSOCIATED(CELLML_MODEL))
        IF(CELLML%MODELS%MODELS(cellml_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          CELLML_MODEL=>CELLML%MODELS%MODELS(cellml_idx)%PTR
        ELSE
          cellml_idx=cellml_idx+1
        ENDIF
      ENDDO
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODEL_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("CELLML_MODEL_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_USER_NUMBER_FIND

  !
  !=================================================================================================================================
  !

  !>Finalises the CellML environments and deallocates all memory.
  SUBROUTINE CELLML_ENVIRONMENTS_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CELLML_ENVIRONMENTS_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    DO WHILE(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS>0)
      CALL CELLML_DESTROY(CELLML_ENVIRONMENTS%ENVIRONMENTS(1)%PTR,ERR,ERROR,*999)
    ENDDO !cellml_idx  

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_ENVIRONMENTS_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_ENVIRONMENTS_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_ENVIRONMENTS_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_ENVIRONMENTS_FINALISE

  !=================================================================================================================================
  !

  !>Initialises the CellML environments.
  SUBROUTINE CELLML_ENVIRONMENTS_INITIALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CELLML_ENVIRONMENTS_INITIALISE",ERR,ERROR,*999)

#ifdef USECELLML

    CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS=0

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_ENVIRONMENTS_INITIALISE")
    RETURN
999 CALL ERRORS("CELLML_ENVIRONMENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_ENVIRONMENTS_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_ENVIRONMENTS_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE CMISS_CELLML
