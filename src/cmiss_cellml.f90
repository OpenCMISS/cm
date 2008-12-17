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

  ! module imports

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
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The (source) field to set up this CellML environment for.
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
  !! - is this where the models are instantiated? (turned into code and compiled)
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
  !! - one field for all models with a component in the field for each variable
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

END MODULE CMISS_CELLML
