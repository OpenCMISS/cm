!> \file
!> $Id: cmiss_cellml.f90 $
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief This module is a openCMISS(cm) buffer module to openCMISS(cellml).
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
!> The Original Code is openCMISS
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

!> This module is a openCMISS(cm) buffer module to openCMISS(cellml).
MODULE CMISS_CELLML

  ! Uses EXAMPLE_PATH from Doxyfile but putting path from examples folder to ensure uniqueness.
  !> \example examples/cellml/src/cellml.f90
  !! A complete example demonstrating and testing the methods defined in the openCMISS(cellml) API.
  !<

  ! module imports
  USE BASE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE TYPES

  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  !PUBLIC public routines go here

CONTAINS

  ! code goes here - just adding it all here to start with, will need to define proper interfaces above.

  !> Set up the CellML environment in the given field.
  !! For a given field, create a CellML environment that will be used to define the value of that field in openCMISS. This will likely simply create and initialise an empty CellML environment object with the specified unique identifier number and associate it with the field? Also set some flag to indicate the CellML environment object is in the process of being created and should not yet be used.
  SUBROUTINE CELLML_CREATE_START(CELLML_USER_NUMBER,FIELD,CELLML,ERR,ERROR,*)
    INTEGER(INTG), INTENT(IN) :: CELLML_USER_NUMBER !<Not sure what the purpose of this is? I guess to give the created CellML environment object a unique identifier? Perhaps if the same identifier is specified then an existing CellML environment would be reused/reinitialised for the given source field...
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The (source) field to set up this CellML environment for. Should be USER_NUMBER instead?
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The newly created CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("CELLML_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_START

  !> Finish creating the CellML environment.
  !! Check the provided CellML environment object and if it all looks good clear the "in progress" flag to indicate the object is now ready for use.
  SUBROUTINE CELLML_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object to check and finialise creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("CELLML_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_FINISH

  !> Prepare the given CellML environment for the creation of individual model objects.
  !! Ensure the provided CellML environment object is suitable and ready for creating all the individual models that will be used with it. Probably simply check the object is valid and empty of all models? or perhaps it also a way to indicate that the user wants to add some more models to it?
  SUBROUTINE CELLML_MODELS_CREATE_START(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object that we want to start creating models with.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("CELLML_MODELS_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_MODELS_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_MODELS_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_CREATE_START

  !> Finalise the given CellML environment after all models have been defined.
  !! This routine is called once all the models for the provided CellML environment have been created.
  !! - is this where the models are instantiated? (turned into code and compiled) - no, need more information before we can generate code.
  !! - perhaps we simply "validate" the models to ensure they are suitable for use before the user gets too far along? or maybe that should be done with the model import routine?
  SUBROUTINE CELLML_MODELS_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object that we want to finalise the creation of new models with.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("CELLML_MODELS_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_MODELS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_MODELS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_CREATE_FINISH

  !> Import the specified CellML model into the given CellML environment object.
  !! Here we load specified CellML models into the CellML environment object. Will be called for each model required for use with the base source field for which this CellML environment was created.
  !! - should do some level of validation when the model is loaded
  !! - see URI notes below...
  SUBROUTINE CELLML_MODEL_IMPORT(MODEL_USER_NUMBER,CELLML,URI,ERR,ERROR,*)
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The unique identifier for this model within the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object into which we want to import the specified model.
    TYPE(VARYING_STRING), INTENT(IN) :: URI !<The (absolute? relative?) URI of the model to import. This should be a fully qualified URI which resolves to either a CellML XML model element or locates CellML simulation metadata (openCMISS CellML related metadata?). In the case of simulation metadata, the model will be completely defined from the metadata. If a bare model is found then the default simulation parameters are used? or the user is prompted for more data? or an error is flagged?
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_MODEL_IMPORT",ERR,ERROR,*999)
    CALL EXITS("CELLML_MODEL_IMPORT")
    RETURN
999 CALL ERRORS("CELLML_MODEL_IMPORT",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_IMPORT")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_IMPORT

  !> Start the creation of the models field for the given CellML environment.
  !! This will create the models field for the given CellML environment object. The models field is used to associate specific models defined for this CellML environment with each of the degrees of freedom for the (source) field for which this CellML environment is defined.
  !! - what to do if models field already exists? exists but has a different user number?
  !! - perhaps a better way to do this is through the mapping of CellML model variables to field components. This might provide more flexibility, including the use of multiple models for different field components at a single DOF. Otherwise, still need some way to map model variables to field components.
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_START(MODEL_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    INTEGER(INTG), INTENT(IN) :: MODEL_FIELD_USER_NUMBER !<The unique identifier for the models field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to create the models field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_CREATE_START

  !> Finish the creation of the models field for the given CellML environment.
  !! This will finalise the creation of the models field for the given CellML environment object.
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to finish creation of the models field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_CREATE_FINISH

  !> Fetch the models field for the given CellML environment.
  !! Check the models field is correctly defined and return it for the user to define. The user will be able to use the returned field object to assign models to DOFs.
  SUBROUTINE CELLML_MODELS_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the models field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the models field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_MODELS_FIELD_GET",ERR,ERROR,*999)
    CALL EXITS("CELLML_MODELS_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_GET

  !> Start the creation of the state field for the given CellML environment.
  !! We want to create a single field which will be used to store the state variables of each of the models in the given CellML environment at each of the DOFs for the environment's base field.
  !! - one field for all models with a component in the field for each variable (perhaps simply use the largest model to set the number of components for all DOFs and some will not be used).
  !! - can a single field have different components at different DOFs (grid points)?
  !! - perhaps need a state field for each model and each field is only defined at the DOFs valid for that model?
  !! - this has to be called after the models field is defined.
  SUBROUTINE CELLML_STATE_FIELD_CREATE_START(STATE_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    INTEGER(INTG), INTENT(IN) :: STATE_FIELD_USER_NUMBER !<The unique identifier for the state field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to create the state field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_STATE_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_CREATE_START

  !> Finialse the creation of the state field for the given CellML environment.
  !! Finish creating the state variable field for the provided CellML environment.
  !! - default field values are set for all components in the field at all DOFs based on the initial_value's from the CellML model
  SUBROUTINE CELLML_STATE_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the creation of the state field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_STATE_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_CREATE_FINISH

  !> Fetch the state field for the given CellML environment.
  !! Check the state field is correctly defined and return it for the user to access - either for getting the current values or to set initial conditions.
  !! - need some way to get from variable URIs in the CellML models to components in the returned state field.
  !! - standard field routines should be used to get/set field component values.
  SUBROUTINE CELLML_STATE_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the state field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the state field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_STATE_FIELD_GET",ERR,ERROR,*999)
    CALL EXITS("CELLML_STATE_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_GET

  !> Find the component ID in the given field for the variable defined by the given URI in the provided CellML environment.
  !! This generic routine will be used to map variable URI's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model URI also needed?
  SUBROUTINE CELLML_FIELD_COMPONENT_GET(CELLML,FIELD,URI,COMPONENT_USER_NUMBER,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The field in which to locate the specific component. No need to know whether this is the state, algebraic, or parameter field?
    TYPE(VARYING_STRING), INTENT(IN) :: URI !<The URI of the model variable which needs to be located in the provided field.
    INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On successful return will be the component user number of the field component for the model variable defined by the given URI.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_FIELD_COMPONENT_GET",ERR,ERROR,*999)
    CALL EXITS("CELLML_FIELD_COMPONENT_GET")
    RETURN
999 CALL ERRORS("CELLML_FIELD_COMPONENT_GET",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_COMPONENT_GET")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_COMPONENT_GET

  !> Create a field used to store intermediate variables of interest.
  !! The intermediate field created for this environment provides access to the value of nominated intermediate variables from the models available in this CellML environment.
  !! - this field will be similar in structure to the state field, but only contains the nominated variables.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_START(INTERMEDIATE_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    INTEGER(INTG), INTENT(IN) :: INTERMEDIATE_FIELD_USER_NUMBER !<The unique identifier for the intermediate field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_START

  !> Finialse the creation of the intermediate field for the given CellML environment.
  !! Finish creating the intermediate variable field for the provided CellML environment.
  !! - check for valid variables being defined?
  !! - maybe delete the field if no valid variables/components exist?
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the creation of the intermediate field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_FINISH

  !> Add a specific variable to this CellML environment's intermediate field.
  !! Nominate a specific variable from a model in this environment for inclusion in the environment's intermediate field. This will ensure the value of this variable is available at each of the environment's (source) field DOF's for which the specified model is valid.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD(CELLML,MODEL_URI,VARIABLE_URI,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to add the nominated variable to the intermediate field.
    TYPE(VARYING_STRING), INTENT(IN) :: MODEL_URI !<The URI of the model (simulation?) in which to look for the variable being nominated.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_URI !<The URI of the variable in the specified model to include in the intermediate field for this CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_ADD",ERR,ERROR,*999)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_ADD",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_ADD")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_ADD

  !> Fetch the intermediate field for the given CellML environment.
  !! Check the intermediate field is correctly defined and return it for the user to access.
  !! - is there a way to ensure the field returned to the calling routine is "read-only"?
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the intermediate field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the intermediate field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_GET",ERR,ERROR,*999)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_GET

  !> Start the parameters definition process.
  !! Initialise the CellML environment ready for the nomination of parameters from this environment's models that will be overridden by field values.
  !! - what to do if the paramters field has already be defined?
  SUBROUTINE CELLML_PARAMETERS_CREATE_START(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we will be defining parameter overrides.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_PARAMETERS_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_PARAMETERS_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_CREATE_START

  !> Finialse the parameters definition process.
  !! Indicates that the user has added all parameter overrides and allows the CellML environment to be further processed.
  SUBROUTINE CELLML_PARAMETERS_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the parameter overrides.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_PARAMETERS_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_PARAMETERS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_CREATE_FINISH

  !> Nominate a specific parameter in CellML environment to override.
  !! Nominate a specific parameter variable from a model in this environment which will have its value overridden by the envronment's parameter field.
  SUBROUTINE CELLML_PARAMETERS_ADD(CELLML,MODEL_URI,VARIABLE_URI,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object.
    TYPE(VARYING_STRING), INTENT(IN) :: MODEL_URI !<The URI of the model (simulation?) in which to look for the variable being nominated.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_URI !<The URI of the variable in the specified model to specify for overriding.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_PARAMETERS_ADD",ERR,ERROR,*999)
    CALL EXITS("CELLML_PARAMETERS_ADD")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_ADD",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_ADD")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_ADD

  !> Create the parameters field for this CellML environment.
  !! Check that parameters have been nominated and create the environment's user accessible parameters field. The user accessible parameters field is used by the user to define parameter overrides and is distinct from the internally used parameters field which is always defined for the environment's (source) field DOF's.
  !! - similar to the models field.
  !! - should initialise all components of the field as constant fields with a value as specified in the corresponding CellML models imported into the environment?
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_START(PARAMETERS_FIELD_USER_NUMBER,CELLML,ERR,ERROR,*)
    INTEGER(INTG), INTENT(IN) :: PARAMETERS_FIELD_USER_NUMBER !<The unique identifier for the parameters field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we will be defining the parameters field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_START",ERR,ERROR,*999)
    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_START

  !> Finialse the parameters field definition.
  !! - here is where we sync the user level parameters field to the internal parameters representation?
  !! - or maybe we simply create the internal field?
  !! - the internal field will contain components for all parameters in all models?
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the parameters field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_FINISH",ERR,ERROR,*999)
    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_FINISH

  !> Fetch the user parameters field for the given CellML environment.
  !! Check the parameters field is correctly defined and return it for the user to access in order to define the desired parameter overrides.
  SUBROUTINE CELLML_PARAMETERS_FIELD_GET(CELLML,FIELD,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the prameters field.
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On successful return will be set to the parameters field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_PARAMETERS_FIELD_GET",ERR,ERROR,*999)
    CALL EXITS("CELLML_PARAMETERS_FIELD_GET")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_GET",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_GET")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_GET

  !> Validate and instantiate the specified CellML environment.
  !! Users should call this routine once they have set up the CellML environment to allow the CellML environment to be validated and instantiated into computable code.
  !! - generate and compile code.
  !! - distribute things amonst processors?
  !! - create the internal fields?
  !! - does this method need to get called or can it be done implicitly as needed by other solution procedures?
  SUBROUTINE CELLML_GENERATE(CELLML,ERR,ERROR,*)
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to generate computable code.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS("CELLML_GENERATE",ERR,ERROR,*999)
    CALL EXITS("CELLML_GENERATE")
    RETURN
999 CALL ERRORS("CELLML_GENERATE",ERR,ERROR)
    CALL EXITS("CELLML_GENERATE")
    RETURN 1
  END SUBROUTINE CELLML_GENERATE

END MODULE CMISS_CELLML
