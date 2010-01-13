!> \file
!> $Id: cmiss_cellml_dummy.f90 542 2009-06-03 17:16:22Z chrispbradley $
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief This module is a dummy OpenCMISS(cm) buffer module to OpenCMISS(cellml).
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

!> This module is a dummy OpenCMISS(cm) buffer module to OpenCMISS(cellml).
MODULE CMISS_CELLML

  !Module imports
  USE BASE_ROUTINES
  USE ISO_VARYING_STRING
  USE INPUT_OUTPUT
  USE KINDS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC CELLML_CREATE_START,CELLML_CREATE_FINISH,CELLML_DESTROY
  
  PUBLIC CELLML_MODELS_CREATE_START,CELLML_MODELS_CREATE_FINISH,CELLML_MODEL_IMPORT

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Destroys the given CellML environment.
  SUBROUTINE CELLML_DESTROY(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local Variables

    CALL ENTERS("CELLML_DESTROY",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_CREATE_START(CELLML_USER_NUMBER,FIELD,CELLML,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CELLML_USER_NUMBER 
    TYPE(FIELD_TYPE), POINTER :: FIELD 
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
 
    CALL ENTERS("CELLML_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
   
    CALL EXITS("CELLML_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finish creating the CellML environment.
  SUBROUTINE CELLML_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL ENTERS("CELLML_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
    CALL EXITS("CELLML_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !> Prepare the given CellML environment for the creation of individual model objects.
  SUBROUTINE CELLML_MODELS_CREATE_START(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object that we want to start creating models with.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL ENTERS("CELLML_MODELS_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
 
    CALL EXITS("CELLML_MODELS_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_MODELS_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finalise the given CellML environment after all models have been defined.
  SUBROUTINE CELLML_MODELS_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL ENTERS("CELLML_MODELS_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
    CALL EXITS("CELLML_MODELS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_MODELS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Import the specified CellML model into the given CellML environment object.
  SUBROUTINE CELLML_MODEL_IMPORT(MODEL_USER_NUMBER,CELLML,URI,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER 
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    CHARACTER(LEN=*) :: URI 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODEL_IMPORT",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

    CALL EXITS("CELLML_MODEL_IMPORT")
    RETURN
999 CALL ERRORS("CELLML_MODEL_IMPORT",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_IMPORT")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_IMPORT

  !
  !=================================================================================================================================
  !

  !>Start the creation of the models field for the given CellML environment.
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_START(MODEL_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: MODEL_FIELD_USER_NUMBER
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to finish creation of the models field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_MODELS_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    TYPE(FIELD_TYPE), POINTER :: FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_GET",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_STATE_FIELD_CREATE_START(STATE_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: STATE_FIELD_USER_NUMBER
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_STATE_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_STATE_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    TYPE(FIELD_TYPE), POINTER :: FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_GET",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
   
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
  SUBROUTINE CELLML_FIELD_COMPONENT_GET(CELLML,FIELD,URI,COMPONENT_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    TYPE(FIELD_TYPE), POINTER :: FIELD 
    TYPE(VARYING_STRING), INTENT(IN) :: URI 
    INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_FIELD_COMPONENT_GET",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    COMPONENT_USER_NUMBER=0
    
    CALL EXITS("CELLML_FIELD_COMPONENT_GET")
    RETURN
999 CALL ERRORS("CELLML_FIELD_COMPONENT_GET",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_COMPONENT_GET")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_COMPONENT_GET

  !
  !=================================================================================================================================
  !

  !>Create a field used to store intermediate variables of interest.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_START(INTERMEDIATE_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: INTERMEDIATE_FIELD_USER_NUMBER
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD(CELLML,MODEL_URI,VARIABLE_URI,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    TYPE(VARYING_STRING), INTENT(IN) :: MODEL_URI
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_URI 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_ADD",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_ADD",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD

  !
  !=================================================================================================================================
  !

  !>Fetch the intermediate field for the given CellML environment.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    TYPE(FIELD_TYPE), POINTER :: FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_GET",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_PARAMETERS_CREATE_START(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_PARAMETERS_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
  
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
  SUBROUTINE CELLML_PARAMETERS_ADD(CELLML,MODEL_URI,VARIABLE_URI,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    TYPE(VARYING_STRING), INTENT(IN) :: MODEL_URI
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_URI 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_ADD",ERR,ERROR,*999)
    
    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

    CALL EXITS("CELLML_PARAMETERS_ADD")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_ADD",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_ADD")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_ADD

  !
  !=================================================================================================================================
  !

  !>Create the parameters field for this CellML environment.
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_START(PARAMETERS_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARAMETERS_FIELD_USER_NUMBER
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_START",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the parameters field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_PARAMETERS_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML 
    TYPE(FIELD_TYPE), POINTER :: FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_GET",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)
    
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
  SUBROUTINE CELLML_GENERATE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_GENERATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("CellML dummy routine. Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

    CALL EXITS("CELLML_GENERATE")
    RETURN
999 CALL ERRORS("CELLML_GENERATE",ERR,ERROR)
    CALL EXITS("CELLML_GENERATE")
    RETURN 1
  END SUBROUTINE CELLML_GENERATE

  !
  !=================================================================================================================================
  !

END MODULE CMISS_CELLML
