!> \file
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
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Chris Bradley
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

  !Module imports
  USE ISO_C_BINDING

  USE BASE_ROUTINES
  
#ifdef USECELLML
  USE CELLML_MODEL_DEFINITION
#endif
  
  ! Moved this usage declaration outside the preprocessor definition,
  ! as this file is included only if CELLML integration is selected.
  ! This fixes problems with the CMAKE FORTRAN parser. Its not detecting
  ! the file (=module) dependency correctly and hence breaks the build
  ! on some platforms.
  USE CMISS_FORTRAN_C

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
  !> CellML model variables being mapped to or from fields will have an initial type of CELLML_UNKNOWN_FIELD. This will be set to
  !> the appropriate type once the model is instantiated and the type can be correctly determined.
  !> \todo Link to appropriate methods for instantiation.
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD = 1 !<CellML models field \see CELLML_FieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_STATE_FIELD = 2 !<CellML state field \see CELLML_FieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_INTERMEDIATE_FIELD = 3 !<CellML intermediate field \see CELLML_FieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_PARAMETERS_FIELD = 4 !<CellML parameters field \see CELLML_FieldTypes,CMISS_CELLML
  !>@}

  !> \addtogroup CELLML_FieldMappingTypes CMISS_CELLML::FieldMappingTypes
  !> \brief CellML field parameter types
  !> \see CMISS_CELLML,OPENCMISS_CellMLFieldTypes
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MAP_TO_FIELD_TYPE = 1 !<A CellML to field mapping type \see CELLML_FieldMappingTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_MAP_FROM_FIELD_TYPE = 2 !<A field to CellML mapping type \see CELLML_FieldMappingTypes,CMISS_CELLML
  !>@}

  !> \addtogroup CELLML_ModelsFieldTypes CMISS_CELLML::ModelsFieldTypes
  !> \brief CellML field parameter types
  !> \see CMISS_CELLML,OPENCMISS_CellMLFieldTypes
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD_NOT_CHECKED = -2!<The CellML environment models field has not been checked. \see CELLML_ModelsFieldTypes,CMISS_CELLML
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD_NOT_CONSTANT =-1 !<The CellML environement models field is not constant. \see CELLML_ModelsFieldTypes,CMISS_CELLML
  !>@}

  !Module types

  !Module variables
 
  !Interfaces

  !> Map a CellML variable type from OpenCMISS(cellml) to a CellML field type \see CELLML_FieldTypes,CMISS_CELLML
  INTERFACE MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE
    MODULE PROCEDURE MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG
  END INTERFACE !MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE

  !> Map a CellML field type to a CellML variable type from OpenCMISS(cellml). \see CELLML_FieldTypes,CMISS_CELLML
  INTERFACE MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE
    MODULE PROCEDURE MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG
  END INTERFACE !MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE

  INTERFACE CELLML_MODEL_IMPORT
    MODULE PROCEDURE CELLML_MODEL_IMPORT_C
    MODULE PROCEDURE CELLML_MODEL_IMPORT_VS
  END INTERFACE !CELLML_MODEL_IMPORT
  
  INTERFACE CELLML_VARIABLE_SET_AS_KNOWN
    MODULE PROCEDURE CELLML_VARIABLE_SET_AS_KNOWN_C
    MODULE PROCEDURE CELLML_VARIABLE_SET_AS_KNOWN_VS
  END INTERFACE !CELLML_VARIABLE_SET_AS_KNOWN

  INTERFACE CELLML_VARIABLE_SET_AS_WANTED
    MODULE PROCEDURE CELLML_VARIABLE_SET_AS_WANTED_C
    MODULE PROCEDURE CELLML_VARIABLE_SET_AS_WANTED_VS
  END INTERFACE !CELLML_VARIABLE_SET_AS_WANTED

  INTERFACE CELLML_CREATE_CELLML_TO_FIELD_MAP
    MODULE PROCEDURE CELLML_CREATE_CELLML_TO_FIELD_MAP_C
    MODULE PROCEDURE CELLML_CREATE_CELLML_TO_FIELD_MAP_VS
  END INTERFACE !CELLML_CREATE_CELLML_TO_FIELD_MAP

  INTERFACE CELLML_CREATE_FIELD_TO_CELLML_MAP
    MODULE PROCEDURE CELLML_CREATE_FIELD_TO_CELLML_MAP_C
    MODULE PROCEDURE CELLML_CREATE_FIELD_TO_CELLML_MAP_VS
  END INTERFACE !CELLML_CREATE_FIELD_TO_CELLML_MAP

  INTERFACE CELLML_FIELD_COMPONENT_GET
    MODULE PROCEDURE CELLML_FIELD_COMPONENT_GET_C
    MODULE PROCEDURE CELLML_FIELD_COMPONENT_GET_VS
  END INTERFACE !CELLML_FIELD_COMPONENT_GET

  PUBLIC MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE,MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE

  PUBLIC CELLML_MODELS_FIELD,CELLML_STATE_FIELD,CELLML_INTERMEDIATE_FIELD,CELLML_PARAMETERS_FIELD

  PUBLIC CELLML_MODELS_FIELD_NOT_CONSTANT

  PUBLIC CELLML_CELLML_TO_FIELD_UPDATE
  
  PUBLIC CELLML_CREATE_START,CELLML_CREATE_FINISH

  PUBLIC CELLML_DESTROY

  PUBLIC CELLML_FIELD_TO_CELLML_UPDATE
  
  PUBLIC CELLML_MODEL_IMPORT

  PUBLIC CELLML_VARIABLE_SET_AS_KNOWN,CELLML_VARIABLE_SET_AS_WANTED

  PUBLIC CELLML_FIELD_MAPS_CREATE_START,CELLML_FIELD_MAPS_CREATE_FINISH

  PUBLIC CellML_FieldModelDofSet
  
  PUBLIC CELLML_CREATE_CELLML_TO_FIELD_MAP,CELLML_CREATE_FIELD_TO_CELLML_MAP

  PUBLIC CELLML_MODELS_FIELD_CREATE_START,CELLML_MODELS_FIELD_CREATE_FINISH,CELLML_MODELS_FIELD_GET

  PUBLIC CELLML_STATE_FIELD_CREATE_START,CELLML_STATE_FIELD_CREATE_FINISH,CELLML_STATE_FIELD_GET

  PUBLIC CELLML_FIELD_COMPONENT_GET

  PUBLIC CELLML_INTERMEDIATE_FIELD_CREATE_FINISH,CELLML_INTERMEDIATE_FIELD_CREATE_START

  PUBLIC CELLML_INTERMEDIATE_FIELD_GET

  PUBLIC CELLML_PARAMETERS_FIELD_CREATE_START,CELLML_PARAMETERS_FIELD_CREATE_FINISH

  PUBLIC CELLML_PARAMETERS_FIELD_GET

  PUBLIC CELLML_GENERATE

  PUBLIC CELLML_USER_NUMBER_FIND

  PUBLIC CELLML_ENVIRONMENTS_FINALISE,CELLML_ENVIRONMENTS_INITIALISE
  
CONTAINS

  !
  !=================================================================================================================================
  !


  !>Updates any mapped fields from the cellml fields
  SUBROUTINE CELLML_CELLML_TO_FIELD_UPDATE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to update
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: derivativeNumber,dofIdx,dofType,dof2ParamIdx,elementNumber,gaussNumber,gridNumber,map_idx,modelIdx, &
      & nodeNumber,versionNumber
    INTEGER(INTG), POINTER :: MODELS_DATA(:)
    REAL(DP) :: dofValue
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: FIELD_MAPS
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: MODEL_MAP
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: MODEL_MAPS
    TYPE(FIELD_TYPE), POINTER :: MODELS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: MODELS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CELLML_CELLML_TO_FIELD_UPDATE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
        FIELD_MAPS=>CELLML%FIELD_MAPS
        IF(ASSOCIATED(FIELD_MAPS)) THEN
          IF(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
            !The CellML environement only uses one model and so we can optimise for this.
            MODEL_MAPS=>FIELD_MAPS%MODEL_MAPS(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX)%PTR
            IF(ASSOCIATED(MODEL_MAPS)) THEN
              IF(DIAGNOSTICS1) THEN
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"CellML to field update:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",CELLML%USER_NUMBER,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  One model index = ",CELLML%MODELS_FIELD% &
                  & ONLY_ONE_MODEL_INDEX,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of model maps = ",MODEL_MAPS% &
                  & NUMBER_OF_FIELDS_MAPPED_TO,ERR,ERROR,*999)
              ENDIF
              !Loop over the number of CellML to field maps
              DO map_idx=1,MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM
                MODEL_MAP=>MODEL_MAPS%FIELDS_MAPPED_FROM(map_idx)%PTR
                IF(ASSOCIATED(MODEL_MAP)) THEN
                  SELECT CASE(MODEL_MAP%CELLML_FIELD_TYPE)
                  CASE(CELLML_MODELS_FIELD)
                    CALL FLAG_ERROR("Cannot map models field.",ERR,ERROR,*999)
                  CASE(CELLML_STATE_FIELD)
                    IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(CELLML%STATE_FIELD%STATE_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                        & MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE,MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER, &
                        & ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("CellML environment state field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(CELLML_INTERMEDIATE_FIELD)
                    IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,MODEL_MAP%FIELD, &
                        & MODEL_MAP%VARIABLE_TYPE,MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("CellML environment intermediate field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(CELLML_PARAMETERS_FIELD)
                    IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,MODEL_MAP%FIELD, &
                        & MODEL_MAP%VARIABLE_TYPE,MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("CellML environment parameters field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The CellML to field model map CellML field type of "// &
                      & TRIM(NUMBER_TO_VSTRING(MODEL_MAP%CELLML_FIELD_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for map index "//TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                      & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                      & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  LOCAL_ERROR="The CellML to field map is not associated for map index "// &
                    & TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                    & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                    & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                IF(DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Map index : ",map_idx,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    CellML field type      = ",MODEL_MAP%CELLML_FIELD_TYPE, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    CellML parameter set   = ",MODEL_MAP%CELLML_PARAMETER_SET, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    CellML variable number = ",MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field user number      = ",MODEL_MAP%FIELD%USER_NUMBER, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field variable type    = ",MODEL_MAP%VARIABLE_TYPE, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field parameter set    = ",MODEL_MAP%FIELD_PARAMETER_SET, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field component number = ",MODEL_MAP%COMPONENT_NUMBER, &
                    & ERR,ERROR,*999)                   
                ENDIF
              ENDDO !map_idx
            ELSE
              LOCAL_ERROR="The CellML field maps models map is not associated for model index "// &
                & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            !More than one model used in the models field.
            MODELS_FIELD=>CELLML%MODELS_FIELD%MODELS_FIELD
            IF(ASSOCIATED(MODELS_FIELD)) THEN
              NULLIFY(MODELS_VARIABLE)
              CALL FIELD_VARIABLE_GET(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,MODELS_VARIABLE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MODELS_DATA,ERR,ERROR,*999)
              IF(DIAGNOSTICS1) THEN
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"CellML to field update:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",CELLML%USER_NUMBER,ERR,ERROR,*999)
              ENDIF
              !Loop over the dofs in the models field
              DO dofIdx=1,MODELS_VARIABLE%TOTAL_NUMBER_OF_DOFS
                modelIdx=MODELS_DATA(dofIdx)
                IF(modelIdx>0) THEN
                  MODEL_MAPS=>FIELD_MAPS%MODEL_MAPS(modelIdx)%PTR
                  IF(ASSOCIATED(MODEL_MAPS)) THEN
                    dofType=MODELS_VARIABLE%DOF_TO_PARAM_MAP%DOF_TYPE(1,dofIdx)
                    dof2ParamIdx=MODELS_VARIABLE%DOF_TO_PARAM_MAP%DOF_TYPE(2,dofIdx)
                    SELECT CASE(dofType)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      !Do nothing
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      elementNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP(1,dof2ParamIdx)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      versionNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(1,dof2ParamIdx)
                      derivativeNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(2,dof2ParamIdx)
                      nodeNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,dof2ParamIdx)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      gridNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%GRID_POINT_DOF2PARAM_MAP(1,dof2ParamIdx)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      gaussNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(1,dof2ParamIdx)
                      elementNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(2,dof2ParamIdx)
                    CASE DEFAULT
                      LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                        & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                    IF(DIAGNOSTICS1) THEN
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  DOF index : ",dofIdx,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Model index : ",modelIdx,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      DOF type = ",dofType,ERR,ERROR,*999)
                      SELECT CASE(dofType)
                      CASE(FIELD_CONSTANT_INTERPOLATION)
                        !Do nothing
                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,ERR,ERROR,*999)
                      CASE(FIELD_NODE_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Version number = ",versionNumber,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Derivative number = ",derivativeNumber, &
                          & ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Node number = ",nodeNumber,ERR,ERROR,*999)
                      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Grid number = ",gridNumber,ERR,ERROR,*999)
                      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss number = ",gaussNumber,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                          & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT                      
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of model maps = ",MODEL_MAPS% &
                        & NUMBER_OF_FIELDS_MAPPED_TO,ERR,ERROR,*999)
                    ENDIF
                    !Loop over the number of CellML to field maps
                    DO map_idx=1,MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM
                      MODEL_MAP=>MODEL_MAPS%FIELDS_MAPPED_FROM(map_idx)%PTR
                      IF(ASSOCIATED(MODEL_MAP)) THEN
                        !Get the CellML field DOF value
                        SELECT CASE(MODEL_MAP%CELLML_FIELD_TYPE)
                        CASE(CELLML_MODELS_FIELD)
                          CALL FLAG_ERROR("Cannot map models field.",ERR,ERROR,*999)
                        CASE(CELLML_STATE_FIELD)
                          IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
                            SELECT CASE(dofType)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_ELEMENT(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,elementNumber,MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue, &
                                & ERR,ERROR,*999)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_NODE(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL Field_ParameterSetGetLocalGaussPoint(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,gaussNumber,elementNumber,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                                & dofValue,ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                                & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ELSE
                            CALL FLAG_ERROR("CellML environment state field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE(CELLML_INTERMEDIATE_FIELD)
                          IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
                            SELECT CASE(dofType)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                                & dofValue,ERR,ERROR,*999)
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_ELEMENT(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_NODE(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL Field_ParameterSetGetLocalGaussPoint(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,gaussNumber,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                                & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ELSE
                            CALL FLAG_ERROR("CellML environment intermediate field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE(CELLML_PARAMETERS_FIELD)
                          IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
                            SELECT CASE(dofType)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                                & dofValue,ERR,ERROR,*999)
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_ELEMENT(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,elementNumber,MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue, &
                                & ERR,ERROR,*999)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_NODE(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL Field_ParameterSetGetLocalGaussPoint(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,gaussNumber,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                                & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT                             
                          ELSE
                            CALL FLAG_ERROR("CellML environment parameters field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The CellML to field model map CellML field type of "// &
                            & TRIM(NUMBER_TO_VSTRING(MODEL_MAP%CELLML_FIELD_TYPE,"*",ERR,ERROR))// &
                            & " is invalid for map index "//TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                            & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                            & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                        !Update the OpenCMISS mapped field DOF value
                        SELECT CASE(dofType)
                        CASE(FIELD_CONSTANT_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,dofValue,ERR,ERROR,*999)
                        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_UPDATE_ELEMENT(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,elementNumber,MODEL_MAP%COMPONENT_NUMBER,dofValue, &
                            & ERR,ERROR,*999)
                        CASE(FIELD_NODE_BASED_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_UPDATE_NODE(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                            & MODEL_MAP%COMPONENT_NUMBER,dofValue,ERR,ERROR,*999)
                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                          CALL Field_ParameterSetUpdateGaussPoint(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,gaussNumber,elementNumber,MODEL_MAP%COMPONENT_NUMBER, &
                            & dofValue,ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                            & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT                     
                      ELSE
                        LOCAL_ERROR="The CellML to field map is not associated for map index "// &
                          & TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                          & TRIM(NUMBER_TO_VSTRING(modelIdx,"*",ERR,ERROR))//" of CellML environment number "// &
                          & TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      IF(DIAGNOSTICS1) THEN
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Map index : ",map_idx,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        CellML field type      = ", &
                          & MODEL_MAP%CELLML_FIELD_TYPE,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        CellML parameter set   = ", &
                          & MODEL_MAP%CELLML_PARAMETER_SET,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        CellML variable number = ", &
                          & MODEL_MAP%CELLML_VARIABLE_NUMBER,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field user number      = ", &
                          & MODEL_MAP%FIELD%USER_NUMBER,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field variable type    = ", &
                          & MODEL_MAP%VARIABLE_TYPE,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field parameter set    = ", &
                          & MODEL_MAP%FIELD_PARAMETER_SET,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field component number = ", &
                          & MODEL_MAP%COMPONENT_NUMBER,ERR,ERROR,*999)                   
                      ENDIF                      
                    ENDDO !map_idx
                  ELSE
                    LOCAL_ERROR="The CellML field maps models map is not associated for model index "// &
                      & TRIM(NUMBER_TO_VSTRING(modelIdx,"*",ERR,ERROR))//" of CellML environment number "// &
                      & TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF !modelIdx>0
              ENDDO !dofIdx              
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MODELS_DATA,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("CellML environment models field models field is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML models field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    END IF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CELLML_TO_FIELD_UPDATE")
    RETURN
999 CALL ERRORS("CELLML_CELLML_TO_FIELD_UPDATE",ERR,ERROR)
    CALL EXITS("CELLML_CELLML_TO_FIELD_UPDATE")
    RETURN 1
  END SUBROUTINE CELLML_CELLML_TO_FIELD_UPDATE

  !
  !=================================================================================================================================
  !

  !>Set up the CellML environment in the given region.
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
    TYPE(CELLML_ENVIRONMENTS_TYPE), POINTER :: CELLML_ENVIRONMENTS
    TYPE(CELLML_PTR_TYPE), ALLOCATABLE :: NEW_CELLML_ENVIRONMENTS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_CELLML)

    CALL ENTERS("CELLML_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(CELLML)) THEN
        CALL FLAG_ERROR("CellML is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(CELLML)
        CALL CELLML_USER_NUMBER_FIND(CELLML_USER_NUMBER,REGION,CELLML,ERR,ERROR,*999)
        IF(ASSOCIATED(CELLML)) THEN
          LOCAL_ERROR="CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML_USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ELSE
          CELLML_ENVIRONMENTS=>REGION%CELLML_ENVIRONMENTS
          IF(ASSOCIATED(CELLML_ENVIRONMENTS)) THEN
            !Allocate new cellml
            NULLIFY(NEW_CELLML)
            CALL CELLML_INITIALISE(NEW_CELLML,ERR,ERROR,*999)
            !Set cellml defaults
            NEW_CELLML%USER_NUMBER = CELLML_USER_NUMBER
            NEW_CELLML%GLOBAL_NUMBER = CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1
            NEW_CELLML%ENVIRONMENTS => CELLML_ENVIRONMENTS
            NEW_CELLML%REGION => REGION
            !Add cellml to the list of cellml environments
            ALLOCATE(NEW_CELLML_ENVIRONMENTS(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new CellML environments.",ERR,ERROR,*999)
            DO cellml_idx=1,CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS
              NEW_CELLML_ENVIRONMENTS(cellml_idx)%PTR=>CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR
            ENDDO !cellml_idx
            NEW_CELLML_ENVIRONMENTS(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1)%PTR=>NEW_CELLML
            CALL MOVE_ALLOC(NEW_CELLML_ENVIRONMENTS,CELLML_ENVIRONMENTS%ENVIRONMENTS)
            CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS=CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS+1
            CELLML=>NEW_CELLML
          ELSE
            CALL FLAG_ERROR("Region CellML environments is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_START")
    RETURN
999 IF(ALLOCATED(NEW_CELLML_ENVIRONMENTS)) DEALLOCATE(NEW_CELLML_ENVIRONMENTS)
    CALL ERRORS("CELLML_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finish creating the CellML environment.
  !>At this point we know all the variables that are known and wanted so can generate the final code.
  !! Check the provided CellML environment object and if it all looks good clear the "in progress" flag to indicate the object is now ready for use.
  SUBROUTINE CELLML_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object to check and finialise creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL
    INTEGER(C_INT) :: ERROR_CODE
    INTEGER(INTG) :: model_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CELLML_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(CELLML%CELLML_FINISHED) THEN
        CALL FLAG_ERROR("CellML environment has already been finished.",ERR,ERROR,*999)
      ELSE
        !Check that we have set up the models
        IF(CELLML%NUMBER_OF_MODELS>0) THEN
          DO model_idx=1,CELLML%NUMBER_OF_MODELS
            !write(*,*) 'model_idx = ',model_idx
            CELLML_MODEL => CELLML%MODELS(model_idx)%PTR
            !CALL CELLML_MODEL_DEFINITION_SET_SAVE_TEMP_FILES(CELLML_MODEL%PTR,1)
            ERROR_CODE = CELLML_MODEL_DEFINITION_INSTANTIATE(CELLML_MODEL%PTR)
            IF(ERROR_CODE /= 0) THEN
              LOCAL_ERROR="Error instantiating CellML model index "//TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
            CELLML_MODEL%NUMBER_OF_STATE = CELLML_MODEL_DEFINITION_GET_N_RATES(CELLML_MODEL%PTR)
            IF(CELLML_MODEL%NUMBER_OF_STATE>CELLML%MAXIMUM_NUMBER_OF_STATE)  &
              & CELLML%MAXIMUM_NUMBER_OF_STATE=CELLML_MODEL%NUMBER_OF_STATE
            CELLML_MODEL%NUMBER_OF_PARAMETERS = CELLML_MODEL_DEFINITION_GET_N_CONSTANTS(CELLML_MODEL%PTR)
            IF(CELLML_MODEL%NUMBER_OF_PARAMETERS>CELLML%MAXIMUM_NUMBER_OF_PARAMETERS)  &
              & CELLML%MAXIMUM_NUMBER_OF_PARAMETERS=CELLML_MODEL%NUMBER_OF_PARAMETERS
            CELLML_MODEL%NUMBER_OF_INTERMEDIATE = CELLML_MODEL_DEFINITION_GET_N_ALGEBRAIC(CELLML_MODEL%PTR)
            IF(CELLML_MODEL%NUMBER_OF_INTERMEDIATE>CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE)  &
              & CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE=CELLML_MODEL%NUMBER_OF_INTERMEDIATE
          ENDDO !model_idx
          CELLML%CELLML_FINISHED = .TRUE.
        ELSE
          CALL FLAG_ERROR("Invalid setup. No models have been imported into the CellML environment.",ERR,ERROR,*999)
        ENDIF
      ENDIF
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
  !> Destroys the given CellML environment.
  SUBROUTINE CELLML_DESTROY(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: cellml_idx,cellml_position
    TYPE(CELLML_ENVIRONMENTS_TYPE), POINTER :: CELLML_ENVIRONMENTS
    TYPE(CELLML_PTR_TYPE), ALLOCATABLE :: NEW_CELLML_ENVIRONMENTS(:)
    
    CALL ENTERS("CELLML_DESTROY",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN

      CELLML_ENVIRONMENTS=>CELLML%ENVIRONMENTS
      IF(ASSOCIATED(CELLML_ENVIRONMENTS)) THEN
        cellml_position=CELLML%GLOBAL_NUMBER
        CALL CELLML_FINALISE(CELLML,ERR,ERROR,*999)
        IF(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS>1) THEN
          ALLOCATE(NEW_CELLML_ENVIRONMENTS(CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocated new CellML environments.",ERR,ERROR,*999)
          DO cellml_idx=1,CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS
            IF(cellml_idx<cellml_position) THEN
              NEW_CELLML_ENVIRONMENTS(cellml_idx)%PTR=>CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR
            ELSE IF(cellml_idx>cellml_position) THEN
              CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR%GLOBAL_NUMBER=CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)% &
                & PTR%GLOBAL_NUMBER-1
              NEW_CELLML_ENVIRONMENTS(cellml_idx-1)%PTR=>CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR              
            ENDIF
          ENDDO !cellml_idx
          CALL MOVE_ALLOC(NEW_CELLML_ENVIRONMENTS,CELLML_ENVIRONMENTS%ENVIRONMENTS)
          CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS=CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS-1
        ELSE
          IF(ALLOCATED(CELLML_ENVIRONMENTS%ENVIRONMENTS)) DEALLOCATE(CELLML_ENVIRONMENTS%ENVIRONMENTS)
          CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS=0
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environments is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML is not associated.",ERR,ERROR,*999)
    END IF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_DESTROY")
    RETURN
999 IF(ALLOCATED(NEW_CELLML_ENVIRONMENTS)) DEALLOCATE(NEW_CELLML_ENVIRONMENTS)
    CALL ERRORS("CELLML_DESTROY",ERR,ERROR)
    CALL EXITS("CELLML_DESTROY")
    RETURN 1
    
  END SUBROUTINE CELLML_DESTROY

  !
  !=================================================================================================================================
  !


  !>Updates any cellml fields from the mapped fields
  SUBROUTINE CELLML_FIELD_TO_CELLML_UPDATE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to update
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: derivativeNumber,dofIdx,dofType,dof2ParamIdx,elementNumber,gaussNumber,gridNumber,map_idx,modelIdx, &
      & nodeNumber,versionNumber
    INTEGER(INTG), POINTER :: MODELS_DATA(:)
    REAL(DP) :: dofValue
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: FIELD_MAPS
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: MODEL_MAP
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: MODEL_MAPS
    TYPE(FIELD_TYPE), POINTER :: MODELS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: MODELS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CELLML_FIELD_TO_CELLML_UPDATE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
        FIELD_MAPS=>CELLML%FIELD_MAPS
        IF(ASSOCIATED(FIELD_MAPS)) THEN
          IF(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
            !The CellML environement only uses one model and so we can optimise for this.
            MODEL_MAPS=>FIELD_MAPS%MODEL_MAPS(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX)%PTR
            IF(ASSOCIATED(MODEL_MAPS)) THEN
              IF(DIAGNOSTICS1) THEN
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Field to CellML update:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",CELLML%USER_NUMBER,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  One model index = ",CELLML%MODELS_FIELD% &
                  & ONLY_ONE_MODEL_INDEX,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of model maps = ",MODEL_MAPS% &
                  & NUMBER_OF_FIELDS_MAPPED_TO,ERR,ERROR,*999)
              ENDIF
             !Loop over the number of field to CellML maps
              DO map_idx=1,MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO
                MODEL_MAP=>MODEL_MAPS%FIELDS_MAPPED_TO(map_idx)%PTR
                IF(ASSOCIATED(MODEL_MAP)) THEN
                  SELECT CASE(MODEL_MAP%CELLML_FIELD_TYPE)
                  CASE(CELLML_MODELS_FIELD)
                    CALL FLAG_ERROR("Cannot map models field.",ERR,ERROR,*999)
                  CASE(CELLML_STATE_FIELD)
                    IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                        & MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,CELLML%STATE_FIELD%STATE_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("CellML environment state field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(CELLML_INTERMEDIATE_FIELD)
                    IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                        & MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("CellML environment intermediate field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(CELLML_PARAMETERS_FIELD)
                    IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
                      CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                        & MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("CellML environment parameters field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The CellML to field model map CellML field type of "// &
                      & TRIM(NUMBER_TO_VSTRING(MODEL_MAP%CELLML_FIELD_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for map index "//TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                      & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                      & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  LOCAL_ERROR="The CellML to field map is not associated for map index "// &
                    & TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                    & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                    & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                IF(DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Map index : ",map_idx,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field user number      = ",MODEL_MAP%FIELD%USER_NUMBER, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field variable type    = ",MODEL_MAP%VARIABLE_TYPE, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field parameter set    = ",MODEL_MAP%FIELD_PARAMETER_SET, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field component number = ",MODEL_MAP%COMPONENT_NUMBER, &
                    & ERR,ERROR,*999)                   
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    CellML field type      = ",MODEL_MAP%CELLML_FIELD_TYPE, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    CellML parameter set   = ",MODEL_MAP%CELLML_PARAMETER_SET, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    CellML variable number = ",MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                    & ERR,ERROR,*999)
                ENDIF
              ENDDO !map_idx
            ELSE
              LOCAL_ERROR="The CellML field maps models map is not associated for model index "// &
                & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            !More than one model is used in the models field
            MODELS_FIELD=>CELLML%MODELS_FIELD%MODELS_FIELD
            IF(ASSOCIATED(MODELS_FIELD)) THEN
              NULLIFY(MODELS_VARIABLE)
              CALL FIELD_VARIABLE_GET(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,MODELS_VARIABLE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MODELS_DATA,ERR,ERROR,*999)
              IF(DIAGNOSTICS1) THEN
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Field to CellML update:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",CELLML%USER_NUMBER,ERR,ERROR,*999)
              ENDIF
              !Loop over the dofs in the models field
              DO dofIdx=1,MODELS_VARIABLE%TOTAL_NUMBER_OF_DOFS
                modelIdx=MODELS_DATA(dofIdx)
                IF(modelIdx>0) THEN
                  MODEL_MAPS=>FIELD_MAPS%MODEL_MAPS(modelIdx)%PTR
                  IF(ASSOCIATED(MODEL_MAPS)) THEN
                    dofType=MODELS_VARIABLE%DOF_TO_PARAM_MAP%DOF_TYPE(1,dofIdx)
                    dof2ParamIdx=MODELS_VARIABLE%DOF_TO_PARAM_MAP%DOF_TYPE(2,dofIdx)
                    SELECT CASE(dofType)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      !Do nothing
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      elementNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP(1,dof2ParamIdx)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      versionNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(1,dof2ParamIdx)
                      derivativeNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(2,dof2ParamIdx)
                      nodeNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,dof2ParamIdx)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      gridNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%GRID_POINT_DOF2PARAM_MAP(1,dof2ParamIdx)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      gaussNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(1,dof2ParamIdx)
                      elementNumber=MODELS_VARIABLE%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(2,dof2ParamIdx)
                    CASE DEFAULT
                      LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                        & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                    IF(DIAGNOSTICS1) THEN
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  DOF index : ",dofIdx,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Model index : ",modelIdx,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      DOF type = ",dofType,ERR,ERROR,*999)
                      SELECT CASE(dofType)
                      CASE(FIELD_CONSTANT_INTERPOLATION)
                        !Do nothing
                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,ERR,ERROR,*999)
                      CASE(FIELD_NODE_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Version number = ",versionNumber,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Derivative number = ",derivativeNumber, &
                          & ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Node number = ",nodeNumber,ERR,ERROR,*999)
                      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Grid number = ",gridNumber,ERR,ERROR,*999)
                      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss number = ",gaussNumber,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                          & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT                      
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of model maps = ",MODEL_MAPS% &
                        & NUMBER_OF_FIELDS_MAPPED_TO,ERR,ERROR,*999)
                    ENDIF
                    !Loop over the number of field to CellML maps
                    DO map_idx=1,MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO
                      MODEL_MAP=>MODEL_MAPS%FIELDS_MAPPED_TO(map_idx)%PTR
                      IF(ASSOCIATED(MODEL_MAP)) THEN
                        !Get the OpenCMISS mapped field DOF value
                        SELECT CASE(dofType)
                        CASE(FIELD_CONSTANT_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_GET_CONSTANT(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,MODEL_MAP%COMPONENT_NUMBER,dofValue,ERR,ERROR,*999)
                        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_GET_ELEMENT(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,elementNumber,MODEL_MAP%COMPONENT_NUMBER,dofValue, &
                            & ERR,ERROR,*999)
                        CASE(FIELD_NODE_BASED_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_GET_NODE(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                            & MODEL_MAP%COMPONENT_NUMBER,dofValue,ERR,ERROR,*999)
                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                          CALL Field_ParameterSetGetLocalGaussPoint(MODEL_MAP%FIELD,MODEL_MAP%VARIABLE_TYPE, &
                            & MODEL_MAP%FIELD_PARAMETER_SET,gaussNumber,elementNumber, &
                            & MODEL_MAP%COMPONENT_NUMBER,dofValue,ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                            & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT                     
                        !Update the CellML field DOF value
                        SELECT CASE(MODEL_MAP%CELLML_FIELD_TYPE)
                        CASE(CELLML_MODELS_FIELD)
                          CALL FLAG_ERROR("Cannot map models field.",ERR,ERROR,*999)
                        CASE(CELLML_STATE_FIELD)
                          IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
                            SELECT CASE(dofType)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_ELEMENT(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,elementNumber,MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue, &
                                & ERR,ERROR,*999)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_NODE(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL Field_ParameterSetUpdateGaussPoint(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,gaussNumber,elementNumber,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                                & dofValue,ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                                & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ELSE
                            CALL FLAG_ERROR("CellML environment state field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE(CELLML_INTERMEDIATE_FIELD)
                          IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
                            SELECT CASE(dofType)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                                & dofValue,ERR,ERROR,*999)
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_ELEMENT(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_NODE(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL Field_ParameterSetUpdateGaussPoint(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,gaussNumber,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                                & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ELSE
                            CALL FLAG_ERROR("CellML environment intermediate field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE(CELLML_PARAMETERS_FIELD)
                          IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
                            SELECT CASE(dofType)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,MODEL_MAP%CELLML_VARIABLE_NUMBER, &
                                & dofValue,ERR,ERROR,*999)
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_ELEMENT(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_UPDATE_NODE(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & MODEL_MAP%CELLML_PARAMETER_SET,versionNumber,derivativeNumber,nodeNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                              CALL Field_ParameterSetUpdateGaussPoint(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,MODEL_MAP%CELLML_PARAMETER_SET,gaussNumber,elementNumber, &
                                & MODEL_MAP%CELLML_VARIABLE_NUMBER,dofValue,ERR,ERROR,*999)
                            CASE DEFAULT
                              LOCAL_ERROR="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",ERR,ERROR))// &
                                & " for local DOF number "//TRIM(NUMBER_TO_VSTRING(dofIdx,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT                             
                          ELSE
                            CALL FLAG_ERROR("CellML environment parameters field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The CellML to field model map CellML field type of "// &
                            & TRIM(NUMBER_TO_VSTRING(MODEL_MAP%CELLML_FIELD_TYPE,"*",ERR,ERROR))// &
                            & " is invalid for map index "//TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                            & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))// &
                            & " of CellML environment number "//TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        LOCAL_ERROR="The CellML to field map is not associated for map index "// &
                          & TRIM(NUMBER_TO_VSTRING(map_idx,"*",ERR,ERROR))//" of model index "// &
                          & TRIM(NUMBER_TO_VSTRING(modelIdx,"*",ERR,ERROR))//" of CellML environment number "// &
                          & TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      IF(DIAGNOSTICS1) THEN
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Map index : ",map_idx,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        CellML field type      = ", &
                          & MODEL_MAP%CELLML_FIELD_TYPE,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        CellML parameter set   = ", &
                          & MODEL_MAP%CELLML_PARAMETER_SET,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        CellML variable number = ", &
                          & MODEL_MAP%CELLML_VARIABLE_NUMBER,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field user number      = ", &
                          & MODEL_MAP%FIELD%USER_NUMBER,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field variable type    = ", &
                          & MODEL_MAP%VARIABLE_TYPE,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field parameter set    = ", &
                          & MODEL_MAP%FIELD_PARAMETER_SET,ERR,ERROR,*999)
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Field component number = ", &
                          & MODEL_MAP%COMPONENT_NUMBER,ERR,ERROR,*999)                   
                      ENDIF                      
                    ENDDO !map_idx
                  ELSE
                    LOCAL_ERROR="The CellML field maps models map is not associated for model index "// &
                      & TRIM(NUMBER_TO_VSTRING(modelIdx,"*",ERR,ERROR))//" of CellML environment number "// &
                      & TRIM(NUMBER_TO_VSTRING(CELLML%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF !modelIdx>0
              ENDDO !dofIdx              
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MODELS_DATA,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("CellML environment models field models field is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment models field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    END IF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FIELD_TO_CELLML_UPDATE")
    RETURN
999 CALL ERRORS("CELLML_FIELD_TO_CELLML_UPDATE",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_TO_CELLML_UPDATE")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_TO_CELLML_UPDATE

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
    INTEGER(INTG) :: model_idx
    
    CALL ENTERS("CELLML_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ALLOCATED(CELLML%MODELS)) THEN
        DO model_idx=1,SIZE(CELLML%MODELS,1)
          CALL CELLML_MODEL_FINALISE(CELLML%MODELS(model_idx)%PTR,ERR,ERROR,*999)
        ENDDO !model_idx
        DEALLOCATE(CELLML%MODELS)
      ENDIF
      CALL CELLML_FIELD_MAPS_FINALISE(CELLML%FIELD_MAPS,ERR,ERROR,*999)
      CALL CELLML_MODELS_FIELD_FINALISE(CELLML%MODELS_FIELD,ERR,ERROR,*999)
      CALL CELLML_STATE_FIELD_FINALISE(CELLML%STATE_FIELD,ERR,ERROR,*999)
      CALL CELLML_INTERMEDIATE_FIELD_FINALISE(CELLML%INTERMEDIATE_FIELD,ERR,ERROR,*999)
      CALL CELLML_PARAMETERS_FIELD_FINALISE(CELLML%PARAMETERS_FIELD,ERR,ERROR,*999)
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
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CALL FLAG_ERROR("CellML environment is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(CELLML,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML environment.",ERR,ERROR,*999)
      CELLML%GLOBAL_NUMBER=0
      CELLML%USER_NUMBER=0
      NULLIFY(CELLML%REGION)
      NULLIFY(CELLML%ENVIRONMENTS)
      CELLML%CELLML_FINISHED=.FALSE.
      CELLML%NUMBER_OF_MODELS=0
      CELLML%MAXIMUM_NUMBER_OF_STATE=0
      CELLML%MAXIMUM_NUMBER_OF_PARAMETERS=0
      CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE=0
      NULLIFY(CELLML%FIELD_MAPS)
      NULLIFY(CELLML%MODELS_FIELD)
      NULLIFY(CELLML%STATE_FIELD)
      NULLIFY(CELLML%INTERMEDIATE_FIELD)
      NULLIFY(CELLML%PARAMETERS_FIELD)
      CELLML%CELLML_GENERATED=.FALSE.
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_INITIALISE")
    RETURN
999 CALL CELLML_FINALISE(CELLML,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Finish creating the field maps for a CellML environment.
  SUBROUTINE CELLML_FIELD_MAPS_CREATE_FINISH(CELLML,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment to finish creating the maps for. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: model_idx
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: CELLML_MODEL_MAPS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("CELLML_FIELD_MAPS_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(CELLML%CELLML_FINISHED) THEN
        IF(ASSOCIATED(CELLML%FIELD_MAPS)) THEN
          IF(CELLML%FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
            CALL FLAG_ERROR("The CellML environment field maps have already been finished.",ERR,ERROR,*999)
          ELSE
            !Check that each model has field mappings
            DO model_idx=1,CELLML%NUMBER_OF_MODELS
              CELLML_MODEL_MAPS=>CELLML%FIELD_MAPS%MODEL_MAPS(model_idx)%PTR
              IF(ASSOCIATED(CELLML_MODEL_MAPS)) THEN
                IF(CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO==0.AND.CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM==0) THEN
                  LOCAL_ERROR="Invalid setup. CellML model index "//TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))// &
                    & " does not have any mappings to or from a field."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="CellML field maps model maps is not associated for model index "// &
                  & TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !model_idx
            CELLML%FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED=.TRUE.
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML is not associated.",ERR,ERROR,*999)
    ENDIF
    
#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FIELD_MAPS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CELLML_FIELD_MAPS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_MAPS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_MAPS_CREATE_FINISH

  !
  !=================================================================================================================================
  !

  !>Start the creation of field maps for a CellML environment.
  SUBROUTINE CELLML_FIELD_MAPS_CREATE_START(CELLML,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment to create the maps for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
  
    CALL ENTERS("CELLML_FIELD_MAPS_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(CELLML%CELLML_FINISHED) THEN
        !Initialise the field maps
        CALL CELLML_FIELD_MAPS_INITIALISE(CELLML,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("CellML environment has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FIELD_MAPS_CREATE_START")
    RETURN
999 CALL ERRORS("CELLML_FIELD_MAPS_CREATE_START",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_MAPS_CREATE_START")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_MAPS_CREATE_START

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML maps and deallocate all memory.
  SUBROUTINE CELLML_FIELD_MAPS_FINALISE(CELLML_FIELD_MAPS,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS !<A pointer to the CellML field maps to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: model_idx
    
    CALL ENTERS("CELLML_FIELD_MAPS_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
      IF(ALLOCATED(CELLML_FIELD_MAPS%MODEL_MAPS)) THEN
        DO model_idx=1,SIZE(CELLML_FIELD_MAPS%MODEL_MAPS,1)
          CALL CELLML_MODEL_MAPS_FINALISE(CELLML_FIELD_MAPS%MODEL_MAPS(model_idx)%PTR,ERR,ERROR,*999)
        ENDDO !model_idx
        DEALLOCATE(CELLML_FIELD_MAPS%MODEL_MAPS)
      ENDIF
      DEALLOCATE(CELLML_FIELD_MAPS)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_FIELD_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_FIELD_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_MAPS_FINALISE

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML field mpas.
  SUBROUTINE CELLML_FIELD_MAPS_INITIALISE(CELLML,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to initialise the field maps for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR,model_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_FIELD_MAPS_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%FIELD_MAPS)) THEN
        CALL FLAG_ERROR("CellML environment field maps is already associated.",ERR,ERROR,*999)
      ELSE
        IF(CELLML%CELLML_FINISHED) THEN
          ALLOCATE(CELLML%FIELD_MAPS,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML field maps.",ERR,ERROR,*999)
          CELLML%FIELD_MAPS%CELLML=>CELLML
          CELLML%FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED=.FALSE.
          NULLIFY(CELLML%FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)
          NULLIFY(CELLML%FIELD_MAPS%SOURCE_FIELD_VARIABLE)
          NULLIFY(CELLML%FIELD_MAPS%SOURCE_FIELD_DOMAIN)
          CELLML%FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE=0
          !CELLML%FIELD_MAPS%NUMBER_OF_SOURCE_DOFS=0
          !CELLML%FIELD_MAPS%TOTAL_NUMBER_OF_SOURCE_DOFS=0
          !CELLML%FIELD_MAPS%GLOBAL_NUMBER_OF_SOURCE_DOFS=0
          ALLOCATE(CELLML%FIELD_MAPS%MODEL_MAPS(CELLML%NUMBER_OF_MODELS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML environment field maps model maps.",ERR,ERROR,*999)
          DO model_idx=1,CELLML%NUMBER_OF_MODELS
            NULLIFY(CELLML%FIELD_MAPS%MODEL_MAPS(model_idx)%PTR)
            CALL CELLML_MODEL_MAPS_INITIALISE(CELLML%FIELD_MAPS%MODEL_MAPS(model_idx)%PTR,ERR,ERROR,*999)
          ENDDO !model_idx
        ELSE
          CALL FLAG_ERROR("CellML environment has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environement is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_FIELD_MAPS_INITIALISE")
    RETURN
999 CALL CELLML_FIELD_MAPS_FINALISE(CELLML%FIELD_MAPS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_FIELD_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_FIELD_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_FIELD_MAPS_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML model and deallocate all memory.
  SUBROUTINE CELLML_MODEL_FINALISE(CELLML_MODEL,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL !<A pointer to the CellML model to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
     
    CALL ENTERS("CELLML_MODEL_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL)) THEN
      IF(C_ASSOCIATED(CELLML_MODEL%PTR)) CALL DESTROY_CELLML_MODEL_DEFINITION(CELLML_MODEL%PTR)
      CELLML_MODEL%MODEL_ID=""
      DEALLOCATE(CELLML_MODEL)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODEL_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_MODEL_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_FINALISE

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML model.
  SUBROUTINE CELLML_MODEL_INITIALISE(CELLML_MODEL,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL !<On return, a pointer to the CellML initialised model. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_MODEL_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL)) THEN
      CALL FLAG_ERROR("CellML model is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(CELLML_MODEL,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML model.",ERR,ERROR,*999)
      NULLIFY(CELLML_MODEL%CELLML)
      CELLML_MODEL%GLOBAL_NUMBER=0
      CELLML_MODEL%MODEL_ID=""
      CELLML_MODEL%PTR = C_NULL_PTR
      CELLML_MODEL%NUMBER_OF_STATE=0
      CELLML_MODEL%NUMBER_OF_INTERMEDIATE=0
      CELLML_MODEL%NUMBER_OF_PARAMETERS=0
   ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_MODEL_INITIALISE")
    RETURN
999 CALL CELLML_MODEL_FINALISE(CELLML_MODEL,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_MODEL_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML model map and deallocate all memory.
  SUBROUTINE CELLML_MODEL_MAP_FINALISE(CELLML_MODEL_MAP,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: CELLML_MODEL_MAP !<A pointer to the CellML model map to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_MODEL_MAP_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL_MAP)) THEN
      CALL ERASE(CELLML_MODEL_MAP%VARIABLE_ID)
      DEALLOCATE(CELLML_MODEL_MAP)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODEL_MAP_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_MODEL_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_MAP_FINALISE

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML model map.
  SUBROUTINE CELLML_MODEL_MAP_INITIALISE(CELLML_MODEL_MAP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: CELLML_MODEL_MAP !<On return, a pointer to the CellML initialised model map field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_MODEL_MAP_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL_MAP)) THEN
      CALL FLAG_ERROR("CellML model map field is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(CELLML_MODEL_MAP,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML model map.",ERR,ERROR,*999)
      CELLML_MODEL_MAP%CELLML_MAP_TYPE = 0
      NULLIFY(CELLML_MODEL_MAP%FIELD)
      CELLML_MODEL_MAP%VARIABLE_TYPE=0
      CELLML_MODEL_MAP%COMPONENT_NUMBER=0
      CELLML_MODEL_MAP%FIELD_PARAMETER_SET=0
      CELLML_MODEL_MAP%VARIABLE_ID=""
      CELLML_MODEL_MAP%CELLML_FIELD_TYPE=0
      CELLML_MODEL_MAP%CELLML_VARIABLE_NUMBER=0
      CELLML_MODEL_MAP%CELLML_PARAMETER_SET=0
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_MODEL_MAP_INITIALISE")
    RETURN
999 CALL CELLML_MODEL_MAP_FINALISE(CELLML_MODEL_MAP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_MODEL_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_MAP_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML model maps and deallocate all memory.
  SUBROUTINE CELLML_MODEL_MAPS_FINALISE(CELLML_MODEL_MAPS,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: CELLML_MODEL_MAPS !<A pointer to the CellML model maps to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: map_idx
    
    CALL ENTERS("CELLML_MODEL_MAPS_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL_MAPS)) THEN
      IF(ALLOCATED(CELLML_MODEL_MAPS%FIELDS_MAPPED_TO)) THEN
        DO map_idx=1,SIZE(CELLML_MODEL_MAPS%FIELDS_MAPPED_TO,1)
          CALL CELLML_MODEL_MAP_FINALISE(CELLML_MODEL_MAPS%FIELDS_MAPPED_TO(map_idx)%PTR,ERR,ERROR,*999)
        ENDDO !map_idx
        DEALLOCATE(CELLML_MODEL_MAPS%FIELDS_MAPPED_TO)
      ENDIF
      IF(ALLOCATED(CELLML_MODEL_MAPS%FIELDS_MAPPED_FROM)) THEN
        DO map_idx=1,SIZE(CELLML_MODEL_MAPS%FIELDS_MAPPED_FROM,1)
          CALL CELLML_MODEL_MAP_FINALISE(CELLML_MODEL_MAPS%FIELDS_MAPPED_FROM(map_idx)%PTR,ERR,ERROR,*999)
        ENDDO !map_idx
        DEALLOCATE(CELLML_MODEL_MAPS%FIELDS_MAPPED_FROM)
      ENDIF
      DEALLOCATE(CELLML_MODEL_MAPS)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODEL_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_MODEL_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_MAPS_FINALISE

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML model maps.
  SUBROUTINE CELLML_MODEL_MAPS_INITIALISE(CELLML_MODEL_MAPS,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: CELLML_MODEL_MAPS !<On return, a pointer to the CellML initialised model maps. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_MODEL_MAPS_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML_MODEL_MAPS)) THEN
      CALL FLAG_ERROR("CellML model maps is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(CELLML_MODEL_MAPS,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML model maps.",ERR,ERROR,*999)
      CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO=0
      CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM=0
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_MODEL_MAPS_INITIALISE")
    RETURN
999 CALL CELLML_MODEL_MAPS_FINALISE(CELLML_MODEL_MAPS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_MODEL_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODEL_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODEL_MAPS_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Import the specified CellML model into the given CellML environment object.
  SUBROUTINE CELLML_MODEL_IMPORT_C(CELLML,URI,MODEL_INDEX,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object into which we want to import the specified model.
    CHARACTER(LEN=*) :: URI !<The (absolute? relative?) URI of the model to import. As per tracker item 2013 comment 8 the URI should now simply point to a CellML document. Can use a relative URL which will be interpreted relative to the CWD of the executed application.
    INTEGER(INTG), INTENT(OUT) :: MODEL_INDEX !<On return, the index for this model within the given CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: cellml_idx
    TYPE(CELLML_MODEL_TYPE), POINTER :: NEW_MODEL
    TYPE(CELLML_MODEL_PTR_TYPE), ALLOCATABLE :: NEW_MODELS(:)
    !INTEGER (C_INT) CODE
    !TYPE(C_PTR) :: CELLML_MODEL
    CHARACTER(256) :: C_URI
    INTEGER(INTG) :: C_URI_L

    CALL ENTERS("CELLML_MODEL_IMPORT_C",ERR,ERROR,*999)

#ifdef USECELLML

    NULLIFY(NEW_MODEL)
    MODEL_INDEX=0

    IF(ASSOCIATED(CELLML)) THEN
      !set up new model object
      !Allocate new CellML model
      NULLIFY(NEW_MODEL)
      CALL CELLML_MODEL_INITIALISE(NEW_MODEL,ERR,ERROR,*999)
      C_URI_L = LEN_TRIM(URI)
      WRITE(C_URI,'(A,A)') URI(1:C_URI_L),C_NULL_CHAR
      NEW_MODEL%PTR = CREATE_CELLML_MODEL_DEFINITION(C_URI)
      IF(C_ASSOCIATED(NEW_MODEL%PTR)) THEN
        NEW_MODEL%GLOBAL_NUMBER = CELLML%NUMBER_OF_MODELS+1
        NEW_MODEL%MODEL_ID=URI(1:C_URI_l)
        !Add model to this CellML environment's list of models
        ALLOCATE(NEW_MODELS(CELLML%NUMBER_OF_MODELS+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new CellML models array.",ERR,ERROR,*999)
        DO cellml_idx=1,CELLML%NUMBER_OF_MODELS
          NEW_MODELS(cellml_idx)%PTR=>CELLML%MODELS(cellml_idx)%PTR
        ENDDO !cellml_idx
        NEW_MODELS(CELLML%NUMBER_OF_MODELS+1)%PTR=>NEW_MODEL
        CALL MOVE_ALLOC(NEW_MODELS,CELLML%MODELS)
        CELLML%NUMBER_OF_MODELS = CELLML%NUMBER_OF_MODELS+1
        MODEL_INDEX=CELLML%NUMBER_OF_MODELS
      ELSE
        CALL FLAG_ERROR("Error instantiating CellML model.",ERR,ERROR,*999)
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
  SUBROUTINE CELLML_MODEL_IMPORT_VS(CELLML,URI,MODEL_INDEX,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object into which we want to import the specified model.
    TYPE(VARYING_STRING), INTENT(IN) :: URI !<The (absolute? relative?) URI of the model to import. As per tracker item 2013 comment 8 the URI should now simply point to a CellML document. Can use a relative URL which will be interpreted relative to the CWD of the executed application.
    INTEGER(INTG), INTENT(OUT) :: MODEL_INDEX !<On return, the index for this model within the given CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODEL_IMPORT_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_MODEL_IMPORT(CELLML,CHAR(URI),MODEL_INDEX,ERR,ERROR,*999)

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

  !>Sets a CellML model variable to be known - i.e., the variable's value will be set by an OpenCMISS field
  SUBROUTINE CELLML_VARIABLE_SET_AS_KNOWN_C(CELLML,MODEL_INDEX,VARIABLE_ID,ERR,ERROR,*)

    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model in which to find the given variable.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The CellML variable to set as known (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: C_NAME
    INTEGER(INTG) :: C_NAME_L
    INTEGER(C_INT) :: ERROR_CODE
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_VARIABLE_SET_AS_KNOWN_C",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(CELLML%CELLML_FINISHED) THEN
        CALL FLAG_ERROR("CellML environment has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(MODEL_INDEX>0.AND.MODEL_INDEX<=CELLML%NUMBER_OF_MODELS) THEN
          CELLML_MODEL=>CELLML%MODELS(MODEL_INDEX)%PTR
          IF(ASSOCIATED(CELLML_MODEL)) THEN
            !All input arguments are ok.
            C_NAME_L = LEN_TRIM(VARIABLE_ID)
            WRITE(C_NAME,'(A,A)') VARIABLE_ID(1:C_NAME_L),C_NULL_CHAR
            ERROR_CODE = CELLML_MODEL_DEFINITION_SET_VARIABLE_AS_KNOWN(CELLML_MODEL%PTR,C_NAME)
            IF(ERROR_CODE /= 0) THEN
              LOCAL_ERROR="The specified variable can not be set as known: "//VARIABLE_ID
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("CellML model is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified model index of "//TRIM(NUMBER_TO_VSTRING(MODEL_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The modex index should be >= 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_VARIABLE_SET_AS_KNOWN_C")
    RETURN
999 CALL ERRORS("CELLML_VARIABLE_SET_AS_KNOWN_C",ERR,ERROR)
    CALL EXITS("CELLML_VARIABLE_SET_AS_KNOWN_C")
    RETURN 1
  END SUBROUTINE CELLML_VARIABLE_SET_AS_KNOWN_C

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be known - i.e., the variable's value will be set by an OpenCMISS field
  SUBROUTINE CELLML_VARIABLE_SET_AS_KNOWN_VS(CELLML,MODEL_USER_NUMBER,VARIABLE_ID,ERR,ERROR,*)

    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The index of the CellML model in which to find the given variable.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_ID !<The CellML variable to set as known (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_VARIABLE_SET_AS_KNOWN_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_VARIABLE_SET_AS_KNOWN(CELLML,MODEL_USER_NUMBER,CHAR(VARIABLE_ID),ERR,ERROR,*999)

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_VARIABLE_SET_AS_KNOWN_VS")
    RETURN
999 CALL ERRORS("CELLML_VARIABLE_SET_AS_KNOWN_VS",ERR,ERROR)
    CALL EXITS("CELLML_VARIABLE_SET_AS_KNOWN_VS")
    RETURN 1
  END SUBROUTINE CELLML_VARIABLE_SET_AS_KNOWN_VS

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be wanted - i.e., the variable's value will used by an OpenCMISS field
  SUBROUTINE CELLML_VARIABLE_SET_AS_WANTED_C(CELLML,MODEL_INDEX,VARIABLE_ID,ERR,ERROR,*)

    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model in which to find the given variable.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The CellML variable to set as wanted (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: C_NAME
    INTEGER(INTG) :: C_NAME_L
    INTEGER(C_INT) :: ERROR_CODE
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_VARIABLE_SET_AS_WANTED_C",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(CELLML%CELLML_FINISHED) THEN
        CALL FLAG_ERROR("CellML environment has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(MODEL_INDEX>0.AND.MODEL_INDEX<=CELLML%NUMBER_OF_MODELS) THEN
          CELLML_MODEL=>CELLML%MODELS(MODEL_INDEX)%PTR
          IF(ASSOCIATED(CELLML_MODEL)) THEN
            !All input arguments are ok.
            C_NAME_L = LEN_TRIM(VARIABLE_ID)
            WRITE(C_NAME,'(A,A)') VARIABLE_ID(1:C_NAME_L),C_NULL_CHAR
            ERROR_CODE = CELLML_MODEL_DEFINITION_SET_VARIABLE_AS_WANTED(CELLML_MODEL%PTR,C_NAME)
            IF(ERROR_CODE .NE. 0) THEN
              LOCAL_ERROR="The specified variable can not be set as wanted: "//VARIABLE_ID
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("CellML model is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified model index of "//TRIM(NUMBER_TO_VSTRING(MODEL_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The model index should be >= 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_VARIABLE_SET_AS_WANTED_C")
    RETURN
999 CALL ERRORS("CELLML_VARIABLE_SET_AS_WANTED_C",ERR,ERROR)
    CALL EXITS("CELLML_VARIABLE_SET_AS_WANTED_C")
    RETURN 1
  END SUBROUTINE CELLML_VARIABLE_SET_AS_WANTED_C

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be wanted - i.e., the variable's value will be used by an OpenCMISS field
  SUBROUTINE CELLML_VARIABLE_SET_AS_WANTED_VS(CELLML,MODEL_USER_NUMBER,VARIABLE_ID,ERR,ERROR,*)

    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The index of the CellML model in which to find the given variable.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_ID !<The CellML variable to set as wanted (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_VARIABLE_SET_AS_WANTED_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_VARIABLE_SET_AS_WANTED(CELLML,MODEL_USER_NUMBER,CHAR(VARIABLE_ID),ERR,ERROR,*999)

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_VARIABLE_SET_AS_WANTED_VS")
    RETURN
999 CALL ERRORS("CELLML_VARIABLE_SET_AS_WANTED_VS",ERR,ERROR)
    CALL EXITS("CELLML_VARIABLE_SET_AS_WANTED_VS")
    RETURN 1
  END SUBROUTINE CELLML_VARIABLE_SET_AS_WANTED_VS

  !
  !=================================================================================================================================
  !

  !>Create a CellML model variable to field variable component map.
  SUBROUTINE CELLML_CREATE_CELLML_TO_FIELD_MAP_C(CELLML,MODEL_INDEX,VARIABLE_ID,CELLML_PARAMETER_SET, &
    & FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,FIELD_PARAMETER_SET,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model to map from.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The ID of the CellML variable in the given model to map from.
    INTEGER(INTG), INTENT(IN) :: CELLML_PARAMETER_SET !<The CellML parameter set to map from.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field to map to.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable to map to.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to map to.
    INTEGER(INTG), INTENT(IN) :: FIELD_PARAMETER_SET !<The field variable parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    CHARACTER(LEN=1,KIND=C_CHAR) :: C_NAME(MAXSTRLEN)
    !INTEGER(INTG) :: C_NAME_L
    INTEGER(C_INT) :: ERROR_C
    INTEGER(INTG) :: CELLML_VARIABLE_TYPE,CELLML_FIELD_TYPE,CELLML_VARIABLE_NUMBER,map_idx
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: NEW_CELLML_MODEL_MAP
    TYPE(CELLML_MODEL_MAP_PTR_TYPE), ALLOCATABLE :: NEW_FIELDS_MAPPED_FROM(:)
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: CELLML_MODEL_MAPS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CELLML_CREATE_CELLML_TO_FIELD_MAP_C",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CALL CMISSF2CString(VARIABLE_ID,C_NAME)
      IF(CELLML%CELLML_FINISHED) THEN
        CELLML_FIELD_MAPS=>CELLML%FIELD_MAPS
        IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
          IF(CELLML_FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
            CALL FLAG_ERROR("CellML field maps have already been finished.",ERR,ERROR,*999)
          ELSE
            NULLIFY(FIELD_VARIABLE)
            CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
            NULLIFY(PARAMETER_SET)
            CALL FIELD_PARAMETER_SET_GET(FIELD,VARIABLE_TYPE,FIELD_PARAMETER_SET,PARAMETER_SET,ERR,ERROR,*999)
            IF(COMPONENT_NUMBER>0.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              IF(MODEL_INDEX>0.AND.MODEL_INDEX<=CELLML%NUMBER_OF_MODELS) THEN
                CELLML_MODEL=>CELLML%MODELS(MODEL_INDEX)%PTR
                IF(ASSOCIATED(CELLML_MODEL)) THEN
                  CELLML_MODEL_MAPS=>CELLML_FIELD_MAPS%MODEL_MAPS(MODEL_INDEX)%PTR
                  IF(ASSOCIATED(CELLML_MODEL_MAPS)) THEN
                    !All input arguments are ok.
                    ! get the type of the variable being mapped
                    ERROR_C = CELLML_MODEL_DEFINITION_GET_VARIABLE_TYPE(CELLML_MODEL%PTR,C_NAME,CELLML_VARIABLE_TYPE)
                    IF(ERROR_C /= 0) THEN
                      LOCAL_ERROR="Failed to get the type of CellML variable: "// &
                      & VARIABLE_ID// &
                      & "; with the error code: "// &
                      & TRIM(NUMBER_TO_VSTRING(ERROR_C,"*",ERR,ERROR))
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    CELLML_FIELD_TYPE=MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE(CELLML_VARIABLE_TYPE,ERR,ERROR)
                    CALL CELLML_FIELD_COMPONENT_GET(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,CELLML_VARIABLE_NUMBER, &
                      & ERR,ERROR,*999)
                    !C_NAME_L = LEN_TRIM(VARIABLE_ID)
                    !WRITE(C_NAME,'(A,A)') C_NAME(1:C_NAME_L),C_NULL_CHAR
                    !CELLML_VARIABLE_NUMBER=CELLML_MODEL_DEFINITION_ADD_MAPPING_TO_FIELD(CELLML_MODEL%PTR,C_NAME)
                    !Now check that the mapped field is consistent with the other mapped fields for the model.
                    IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)) THEN
                      IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,FIELD%GEOMETRIC_FIELD)) THEN
                        LOCAL_ERROR="The geometric field for field user number "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the geometric field for other field variable components mapped" // &
                          & " in the CellML environment."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN, &
                        & FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN)) THEN
                        LOCAL_ERROR="The domain for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field user number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the domain for other field variable components mapped in the CellML environment."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      IF(CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE/= &
                        & FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE) THEN
                        LOCAL_ERROR="The interpolation type of "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE, &
                          & "*",ERR,ERROR))//" for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_vSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field user number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the interpolation type of "// &
                          & TRIM(NUMBER_TO_VSTRING(CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                          & " used in other field variable components mapped in the CellML environment."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD=>FIELD%GEOMETRIC_FIELD
                      CELLML_FIELD_MAPS%SOURCE_FIELD_VARIABLE=>FIELD_VARIABLE
                      CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN=>FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN
                      CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & INTERPOLATION_TYPE
                    ENDIF
                    !Everything is OK so create the model map field.
                    ! get the type of the variable being mapped
                    ERROR_C = CELLML_MODEL_DEFINITION_GET_VARIABLE_TYPE(CELLML_MODEL%PTR,C_NAME,CELLML_VARIABLE_TYPE)
                    IF(ERROR_C /= 0) THEN
                      LOCAL_ERROR="Failed to get the type of CellML variable: "// &
                      & VARIABLE_ID// &
                      & "; with the error code: "// &
                      & TRIM(NUMBER_TO_VSTRING(ERROR_C,"*",ERR,ERROR))
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    CELLML_FIELD_TYPE=MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE(CELLML_VARIABLE_TYPE,ERR,ERROR)
                    CALL CELLML_FIELD_COMPONENT_GET(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,CELLML_VARIABLE_NUMBER, &
                      & ERR,ERROR,*999)
                    NULLIFY(NEW_CELLML_MODEL_MAP)
                    CALL CELLML_MODEL_MAP_INITIALISE(NEW_CELLML_MODEL_MAP,ERR,ERROR,*999)
                    NEW_CELLML_MODEL_MAP%CELLML_MAP_TYPE=CELLML_MAP_FROM_FIELD_TYPE
                    NEW_CELLML_MODEL_MAP%FIELD=>FIELD
                    NEW_CELLML_MODEL_MAP%VARIABLE_TYPE=VARIABLE_TYPE
                    NEW_CELLML_MODEL_MAP%COMPONENT_NUMBER=COMPONENT_NUMBER
                    NEW_CELLML_MODEL_MAP%FIELD_PARAMETER_SET=FIELD_PARAMETER_SET
                    NEW_CELLML_MODEL_MAP%VARIABLE_ID=VARIABLE_ID
                    NEW_CELLML_MODEL_MAP%CELLML_FIELD_TYPE=CELLML_FIELD_TYPE
                    NEW_CELLML_MODEL_MAP%CELLML_VARIABLE_NUMBER=CELLML_VARIABLE_NUMBER
                    NEW_CELLML_MODEL_MAP%CELLML_PARAMETER_SET=CELLML_PARAMETER_SET
                    !Put this model map field into the list of to field maps
                    ALLOCATE(NEW_FIELDS_MAPPED_FROM(CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM+1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new fields mapped from.",ERR,ERROR,*999)
                    DO map_idx=1,CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM
                      NEW_FIELDS_MAPPED_FROM(map_idx)%PTR=>CELLML_MODEL_MAPS%FIELDS_MAPPED_FROM(map_idx)%PTR
                    ENDDO !map_idx
                    NEW_FIELDS_MAPPED_FROM(CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM+1)%PTR=>NEW_CELLML_MODEL_MAP
                    CALL MOVE_ALLOC(NEW_FIELDS_MAPPED_FROM,CELLML_MODEL_MAPS%FIELDS_MAPPED_FROM)
                    CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM=CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_FROM+1          
                    
                    IF(DIAGNOSTICS1) THEN
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"CellML model variable -> field map:",ERR,ERROR,*999)
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," CellML model :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML User number      = ",CELLML%USER_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Model index             = ",MODEL_INDEX,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Variable ID             = ",VARIABLE_ID,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML field type       = ",CELLML_FIELD_TYPE, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML variable number  = ",CELLML_VARIABLE_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML parameter set    = ",CELLML_PARAMETER_SET, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," Field :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   User number             = ",FIELD%USER_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Variable type           = ",VARIABLE_TYPE, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Component number        = ",COMPONENT_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Parameter set           = ",FIELD_PARAMETER_SET, &
                        & ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("CellML field maps model maps is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("CellML model is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The specified model index of "//TRIM(NUMBER_TO_VSTRING(MODEL_INDEX,"*",ERR,ERROR))// &
                  & " is invalid. The modex index should be >= 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_CELLML_TO_FIELD_MAP_C")
    RETURN
999 CALL ERRORS("CELLML_CREATE_CELLML_TO_FIELD_MAP_C",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_CELLML_TO_FIELD_MAP_C")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_CELLML_TO_FIELD_MAP_C

  !
  !=================================================================================================================================
  !

  !>Create a CellML model variable to field variable component map.
  SUBROUTINE CELLML_CREATE_CELLML_TO_FIELD_MAP_VS(CELLML,MODEL_USER_NUMBER,VARIABLE_ID,CELLML_PARAMETER_SET, &
    & FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,FIELD_PARAMETER_SET,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the CellML model to map from.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_ID !<The ID of the CellML variable in the given model to map from.
    INTEGER(INTG), INTENT(IN) :: CELLML_PARAMETER_SET !<The CellML parameter set to map from.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field to map to.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable to map to.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to map to.
    INTEGER(INTG), INTENT(IN) :: FIELD_PARAMETER_SET !<The field variable parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_CREATE_CELLML_TO_FIELD_MAP_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_CREATE_CELLML_TO_FIELD_MAP(CELLML,MODEL_USER_NUMBER,CHAR(VARIABLE_ID),CELLML_PARAMETER_SET, &
      & FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,FIELD_PARAMETER_SET,ERR,ERROR,*999)

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_CREATE_CELLML_TO_FIELD_MAP_VS")
    RETURN
999 CALL ERRORS("CELLML_CREATE_CELLML_TO_FIELD_MAP_VS",ERR,ERROR)
    CALL EXITS("CELLML_CREATE_CELLML_TO_FIELD_MAP_VS")
    RETURN 1
  END SUBROUTINE CELLML_CREATE_CELLML_TO_FIELD_MAP_VS

  !
  !=================================================================================================================================
  !

  !>Create a field variable component to CellML model variable map.
  SUBROUTINE CELLML_CREATE_FIELD_TO_CELLML_MAP_C(CELLML,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,FIELD_PARAMETER_SET, &
    & MODEL_INDEX,VARIABLE_ID,CELLML_PARAMETER_SET,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field to map from.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to map from.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to map from.
    INTEGER(INTG), INTENT(IN) :: FIELD_PARAMETER_SET !<The field variable parameter set to map from.
    INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model to map to.
    CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The ID of the CellML variable in the given model to map to.
    INTEGER(INTG), INTENT(IN) :: CELLML_PARAMETER_SET !<The CellML parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: CELLML_FIELD_TYPE,CELLML_VARIABLE_NUMBER,map_idx
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS
    TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL
    TYPE(CELLML_MODEL_MAP_TYPE), POINTER :: NEW_CELLML_MODEL_MAP
    TYPE(CELLML_MODEL_MAP_PTR_TYPE), ALLOCATABLE :: NEW_FIELDS_MAPPED_TO(:)
    TYPE(CELLML_MODEL_MAPS_TYPE), POINTER :: CELLML_MODEL_MAPS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    CHARACTER(LEN=1,KIND=C_CHAR) :: C_NAME(MAXSTRLEN)
    !INTEGER(INTG) :: C_NAME_L
    INTEGER(C_INT) :: ERROR_C
    INTEGER(INTG) :: CELLML_VARIABLE_TYPE
    
    CALL ENTERS("CELLML_CREATE_FIELD_TO_CELLML_MAP_C",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CALL CMISSF2CString(VARIABLE_ID,C_NAME)
      IF(CELLML%CELLML_FINISHED) THEN
        CELLML_FIELD_MAPS=>CELLML%FIELD_MAPS
        IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
          IF(CELLML_FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
            CALL FLAG_ERROR("CellML field maps have already been finished.",ERR,ERROR,*999)
          ELSE
            NULLIFY(FIELD_VARIABLE)
            CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
            NULLIFY(PARAMETER_SET)
            CALL FIELD_PARAMETER_SET_GET(FIELD,VARIABLE_TYPE,FIELD_PARAMETER_SET,PARAMETER_SET,ERR,ERROR,*999)
            IF(COMPONENT_NUMBER>0.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              IF(MODEL_INDEX>0.AND.MODEL_INDEX<=CELLML%NUMBER_OF_MODELS) THEN
                CELLML_MODEL=>CELLML%MODELS(MODEL_INDEX)%PTR
                IF(ASSOCIATED(CELLML_MODEL)) THEN
                  CELLML_MODEL_MAPS=>CELLML_FIELD_MAPS%MODEL_MAPS(MODEL_INDEX)%PTR
                  IF(ASSOCIATED(CELLML_MODEL_MAPS)) THEN
                    !All input arguments are ok.
                    ! get the type of the variable being mapped
                    ERROR_C = CELLML_MODEL_DEFINITION_GET_VARIABLE_TYPE(CELLML_MODEL%PTR,C_NAME,CELLML_VARIABLE_TYPE)
                    IF(ERROR_C /= 0) THEN
                      LOCAL_ERROR="Failed to get the type of CellML variable: "// &
                      & VARIABLE_ID// &
                      & "; with the error code: "// &
                      & TRIM(NUMBER_TO_VSTRING(ERROR_C,"*",ERR,ERROR))
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    CELLML_FIELD_TYPE=MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE(CELLML_VARIABLE_TYPE,ERR,ERROR)
                    CALL CELLML_FIELD_COMPONENT_GET(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,CELLML_VARIABLE_NUMBER, &
                      & ERR,ERROR,*999)
                    !C_NAME_L = LEN_TRIM(VARIABLE_ID)
                    !WRITE(C_NAME,'(A,A)') C_NAME(1:C_NAME_L),C_NULL_CHAR
                    !CELLML_VARIABLE_NUMBER=CELLML_MODEL_DEFINITION_ADD_MAPPING_TO_FIELD(CELLML_MODEL%PTR,C_NAME)
                    !Now check that the mapped field is consistent with the other mapped fields for the model.
                    IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)) THEN
                      IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,FIELD%GEOMETRIC_FIELD)) THEN
                        LOCAL_ERROR="The geometric field for field user number "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the geometric field for other field variable components mapped in the" // &
                          & " CellML environment."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN, &
                        & FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN)) THEN
                        LOCAL_ERROR="The domain for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field user number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the domain for other field variable components mapped in the CellML environment."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      IF(CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE/= &
                        & FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE) THEN
                        LOCAL_ERROR="The interpolation type of "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE, &
                          & "*",ERR,ERROR))//" for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field user number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the interpolation type of "// &
                          & TRIM(NUMBER_TO_VSTRING(CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                          & " used in other field variable components mapped in the CellML environment."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN=>FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN
                      CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & INTERPOLATION_TYPE
                    ENDIF
                    !Everything is OK so create the model map field.
                    NULLIFY(NEW_CELLML_MODEL_MAP)
                    CALL CELLML_MODEL_MAP_INITIALISE(NEW_CELLML_MODEL_MAP,ERR,ERROR,*999)
                    NEW_CELLML_MODEL_MAP%CELLML_MAP_TYPE=CELLML_MAP_TO_FIELD_TYPE
                    NEW_CELLML_MODEL_MAP%FIELD=>FIELD
                    NEW_CELLML_MODEL_MAP%VARIABLE_TYPE=VARIABLE_TYPE
                    NEW_CELLML_MODEL_MAP%COMPONENT_NUMBER=COMPONENT_NUMBER
                    NEW_CELLML_MODEL_MAP%FIELD_PARAMETER_SET=FIELD_PARAMETER_SET
                    NEW_CELLML_MODEL_MAP%VARIABLE_ID=VARIABLE_ID
                    NEW_CELLML_MODEL_MAP%CELLML_FIELD_TYPE=CELLML_FIELD_TYPE
                    NEW_CELLML_MODEL_MAP%CELLML_VARIABLE_NUMBER=CELLML_VARIABLE_NUMBER
                    NEW_CELLML_MODEL_MAP%CELLML_PARAMETER_SET=CELLML_PARAMETER_SET
                    !Put this model map field into the list of to field maps
                    ALLOCATE(NEW_FIELDS_MAPPED_TO(CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO+1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new fields mapped to.",ERR,ERROR,*999)
                    DO map_idx=1,CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO
                      NEW_FIELDS_MAPPED_TO(map_idx)%PTR=>CELLML_MODEL_MAPS%FIELDS_MAPPED_TO(map_idx)%PTR
                    ENDDO !map_idx
                    NEW_FIELDS_MAPPED_TO(CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO+1)%PTR=>NEW_CELLML_MODEL_MAP
                    CALL MOVE_ALLOC(NEW_FIELDS_MAPPED_TO,CELLML_MODEL_MAPS%FIELDS_MAPPED_TO)
                    CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO=CELLML_MODEL_MAPS%NUMBER_OF_FIELDS_MAPPED_TO+1
                    
                    IF(DIAGNOSTICS1) THEN
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"CellML model field -> CellML map:",ERR,ERROR,*999)
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," Field :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   User number             = ",FIELD%USER_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Variable type           = ",VARIABLE_TYPE,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Component number        = ",COMPONENT_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Parameter set           = ",FIELD_PARAMETER_SET, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," CellML model :",ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML User number      = ",CELLML%USER_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Model index             = ",MODEL_INDEX,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Variable ID             = ",VARIABLE_ID,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML field type       = ",CELLML_FIELD_TYPE, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML variable number  = ",CELLML_VARIABLE_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   CellML parameter set    = ",CELLML_PARAMETER_SET, &
                        & ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("CellML field maps model maps is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("CellML model is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The specified model index of "//TRIM(NUMBER_TO_VSTRING(MODEL_INDEX,"*",ERR,ERROR))// &
                  & " is invalid. The modex index should be >= 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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
  SUBROUTINE CELLML_CREATE_FIELD_TO_CELLML_MAP_VS(CELLML,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,FIELD_PARAMETER_SET, &
    & MODEL_USER_NUMBER,VARIABLE_ID,CELLML_PARAMETER_SET,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object in which to create the map.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field to map from.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to map from.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to map from.
    INTEGER(INTG), INTENT(IN) :: FIELD_PARAMETER_SET !<The field variable parameter set to map from.
    INTEGER(INTG), INTENT(IN) :: MODEL_USER_NUMBER !<The user number of the CellML model to map to.
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_ID !<The ID of the CellML variable in the given model to map to.
    INTEGER(INTG), INTENT(IN) :: CELLML_PARAMETER_SET !<The CellML parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_CREATE_FIELD_TO_CELLML_MAP_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_CREATE_FIELD_TO_CELLML_MAP(CELLML,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,FIELD_PARAMETER_SET, &
      & MODEL_USER_NUMBER,CHAR(VARIABLE_ID),CELLML_PARAMETER_SET,ERR,ERROR,*999)
    
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

  !>Set the dof in a field specified by a model DOF and component to a value.
  SUBROUTINE CellML_FieldModelDofSet(modelVariable,modelDofIdx,field,variableType,parameterSetIdx,componentIdx, &
     value,err,error,*)
    
    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: modelVariable !<A pointer to the model field variable
    INTEGER(INTG), INTENT(IN) :: modelDofIdx !<The model dof to set.
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to set the value for
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the value for
    INTEGER(INTG), INTENT(IN) :: parameterSetIdx !<The parameter set index of the field variable to set.
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to set.
    REAL(DP), INTENT(IN) :: value !<The value to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: derivativeNumber,dofParamIdx,dofType,elementNumber,gaussNumber,gridNumber,nodeNumber,versionNumber
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("CellML_FieldModelDofSet",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(modelVariable)) THEN
      IF(modelDofIdx>0.AND.modelDofIdx<=modelVariable%TOTAL_NUMBER_OF_DOFS) THEN
        dofType=modelVariable%DOF_TO_PARAM_MAP%DOF_TYPE(1,modelDofIdx)
        dofParamIdx=modelVariable%DOF_TO_PARAM_MAP%DOF_TYPE(2,modelDofIdx)
        SELECT CASE(dofType)
        CASE(FIELD_CONSTANT_DOF_TYPE)
          CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(field,variableType,parameterSetIdx,componentIdx,VALUE,err,error,*999)
        CASE(FIELD_ELEMENT_DOF_TYPE)
          elementNumber=modelVariable%DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP(1,dofParamIdx)
          CALL FIELD_PARAMETER_SET_UPDATE_ELEMENT(field,variableType,parameterSetIdx,elementNumber,componentIdx,value, &
            & err,error,*999)
        CASE(FIELD_NODE_DOF_TYPE)
          versionNumber=modelVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(1,dofParamIdx)
          derivativeNumber=modelVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(2,dofParamIdx)
          nodeNumber=modelVariable%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,dofParamIdx)
          CALL FIELD_PARAMETER_SET_UPDATE_NODE(field,variableType,parameterSetIdx,versionNumber,derivativeNumber,NodeNumber, &
            & componentIdx,VALUE,err,error,*999)
        CASE(FIELD_GRID_POINT_DOF_TYPE)
          gridNumber=modelVariable%DOF_TO_PARAM_MAP%GRID_POINT_DOF2PARAM_MAP(1,dofParamIdx)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_DOF_TYPE)
          gaussNumber=modelVariable%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(1,dofParamIdx)
          elementNumber=modelVariable%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(2,dofParamIdx)
          CALL Field_ParameterSetUpdateGaussPoint(field,variableType,parameterSetIdx,gaussNumber,elementNumber,componentIdx, &
            & VALUE,err,error,*999)
        CASE DEFAULT
          localError="The DOF type of "//TRIM(NUMBER_TO_VSTRING(dofType,"*",err,error))// &
            & " for DOF number "//TRIM(NUMBER_TO_VSTRING(modelDofIdx,"*",err,error))//" is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
       END SELECT
      ELSE
        localError="The model DOF index of "//TRIM(NUMBER_TO_VSTRING(modelDofIdx,"*",err,error))// &
          & " is invalid. The DOF index needs to be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(modelVariable%TOTAL_NUMBER_OF_DOFS,"*",err,error))//"."
        CALL FLAG_ERROR(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Model variable is not asssociated.",err,error,*999)
    ENDIF
    
#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",err,error,*999)

#endif

    CALL EXITS("CellML_FieldModelDofSet")
    RETURN
999 CALL ERRORS("CellML_FieldModelDofSet",err,error)
    CALL EXITS("CellML_FieldModelDofSet")
    RETURN 1
    
  END SUBROUTINE CellML_FieldModelDofSet

  !
  !=================================================================================================================================
  !

  !>Checks a CellML environment models field for correctness.
  SUBROUTINE CELLML_MODELS_FIELD_CHECK(MODELS_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_MODELS_FIELD_TYPE), POINTER :: MODELS_FIELD !<A pointer to the CellML environment models field to check.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: model_idx,source_dof_idx,first_dof_idx
    INTEGER(INTG), POINTER :: MODELS_DATA(:)
    TYPE(CELLML_TYPE), POINTER :: CELLML
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: MODELS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CELLML_MODELS_FIELD_CHECK",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(MODELS_FIELD)) THEN
      IF(MODELS_FIELD%MODELS_FIELD_FINISHED) THEN
        IF(MODELS_FIELD%ONLY_ONE_MODEL_INDEX==CELLML_MODELS_FIELD_NOT_CHECKED) THEN
          CELLML=>MODELS_FIELD%CELLML
          IF(ASSOCIATED(CELLML)) THEN
            !Models field has not been checked before.
            NULLIFY(MODELS_VARIABLE)
            CALL FIELD_VARIABLE_GET(MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,MODELS_VARIABLE,ERR,ERROR,*999)
            IF(MODELS_VARIABLE%NUMBER_OF_DOFS>0) THEN
              CALL FIELD_PARAMETER_SET_DATA_GET(MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MODELS_DATA,ERR,ERROR,*999)
              !check for the first non-zero model index
              source_dof_idx=0
              first_dof_idx=1
              DO source_dof_idx=1,MODELS_VARIABLE%TOTAL_NUMBER_OF_DOFS
                model_idx=MODELS_DATA(source_dof_idx)
                IF(model_idx>=0) THEN
                  MODELS_FIELD%ONLY_ONE_MODEL_INDEX=model_idx
                  first_dof_idx=source_dof_idx
                  IF(model_idx>0) THEN
                    EXIT
                  ENDIF
                ELSE
                  LOCAL_ERROR="The model index of "//TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))// &
                   & " is invalid for source DOF 1. The model index must be >= 0 and <= "// &
                   & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR("The models field has not been set for DOF 1.",ERR,ERROR,*999)
                ENDIF
              ENDDO
              IF(model_idx>=0.AND.model_idx<=CELLML%NUMBER_OF_MODELS) THEN
                DO source_dof_idx=(first_dof_idx+1),MODELS_VARIABLE%TOTAL_NUMBER_OF_DOFS
                  model_idx=MODELS_DATA(source_dof_idx)
                  IF(model_idx>=0.AND.model_idx<=CELLML%NUMBER_OF_MODELS) THEN
                    IF(model_idx/=MODELS_FIELD%ONLY_ONE_MODEL_INDEX.AND.model_idx/=0) THEN
                      MODELS_FIELD%ONLY_ONE_MODEL_INDEX=CELLML_MODELS_FIELD_NOT_CONSTANT
                      EXIT
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The model index of "//TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))// &
                      & " is invalid for source DOF "//TRIM(NUMBER_TO_VSTRING(source_dof_idx,"*",ERR,ERROR))// &
                      & ". The model index must be >= 0 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !source_dof_idx
                IF(MODELS_FIELD%ONLY_ONE_MODEL_INDEX==0) &
                  & CALL FLAG_ERROR("Models field does not have any models set.",ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="The model index of "//TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))// &
                  & " is invalid for source DOF 1. The model index must be >= 0 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR("The models field has not been set for DOF 1.",ERR,ERROR,*999)
              ENDIF
!!TODO: Do we need to make sure it is the same model number on different ranks? The only one model optimisation is to ensure
!!that we don't have to reference the models field inside dof loops on the rank??? 
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MODELS_DATA,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("CellML models field variable does not have any DOFs.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The models field CellML environment is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Models field has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Models field is not associated.",ERR,ERROR,*999)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODELS_FIELD_CHECK")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_CHECK",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_CHECK")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_CHECK

  !
  !=================================================================================================================================
  !

  !>Start the creation of the models field for the given CellML environment.
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_START(MODEL_FIELD_USER_NUMBER,CELLML,MODELS_FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: MODEL_FIELD_USER_NUMBER !<The unique identifier for the models field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to create the models field.
    TYPE(FIELD_TYPE), POINTER :: MODELS_FIELD !<If associated on entry, a pointer to the user created models field which has the same user number as the specified models field user number. If not associated on entry, on exit, a pointer to the created models field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,MODELS_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CELLML_FIELD_MAPS=>CELLML%FIELD_MAPS
      IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
        IF(CELLML_FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
          IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
            CALL FLAG_ERROR("The CellML environment models field is already associated.",ERR,ERROR,*999)
          ELSE
            REGION=>CELLML%REGION
            IF(ASSOCIATED(REGION)) THEN
              IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)) THEN
                IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN)) THEN
                  IF(ASSOCIATED(MODELS_FIELD)) THEN
                    !Check the field has been finished
                    IF(MODELS_FIELD%FIELD_FINISHED) THEN
                      !Check the user numbers match
                      IF(MODEL_FIELD_USER_NUMBER/=MODELS_FIELD%USER_NUMBER) THEN
                        LOCAL_ERROR="The specified models field user number of "// &
                          & TRIM(NUMBER_TO_VSTRING(MODEL_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the user number of the specified models field of "// &
                          & TRIM(NUMBER_TO_VSTRING(MODELS_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      MODELS_FIELD_REGION=>MODELS_FIELD%REGION
                      IF(ASSOCIATED(MODELS_FIELD_REGION)) THEN
                        !Check the field is defined on the same region as the CellML region
                        IF(MODELS_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                          LOCAL_ERROR="Invalid region setup. The specified models field has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(MODELS_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & " and the specified CellML environment has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                        !Check the specified models field has the same geometric field as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,MODELS_FIELD%GEOMETRIC_FIELD)) THEN
                          CALL FLAG_ERROR("The specified models field does not have the same geometric field as the "// &
                            & "geometric field for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                        !Check the specified models field has the same decomposition as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,MODELS_FIELD%DECOMPOSITION)) THEN
                          CALL FLAG_ERROR("The specified models field does not have the same decomposition as the source "// &
                            & "domain decomposition for the specified CellML environment.",ERR,ERROR,*999)
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
                  CALL CELLML_MODELS_FIELD_INITIALISE(CELLML,ERR,ERROR,*999)
                  IF(ASSOCIATED(MODELS_FIELD)) THEN
                    !Now check the supplied field.
                    CALL FIELD_DATA_TYPE_CHECK(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,ERR,ERROR,*999)
                    CALL FIELD_TYPE_CHECK(MODELS_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_CHECK(MODELS_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_CHECK(MODELS_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_CHECK(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                  ELSE
                    CELLML%MODELS_FIELD%MODELS_FIELD_AUTO_CREATED=.TRUE.
                    !Create the CellML environment models field
                    CALL FIELD_CREATE_START(MODEL_FIELD_USER_NUMBER,REGION,CELLML%MODELS_FIELD%MODELS_FIELD,ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE, &
                      & ERR,ERROR,*999)
                    CALL FIELD_LABEL_SET(CELLML%MODELS_FIELD%MODELS_FIELD,"CellMLModelsField",ERR,ERROR,*999)
                    CALL FIELD_TYPE_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,ERR,ERROR,*999)
                    CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,CELLML_FIELD_MAPS% &
                      & SOURCE_GEOMETRIC_FIELD,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_LABEL_SET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,"ModelMap",ERR,ERROR,*999)
                    CALL FIELD_DOF_ORDER_TYPE_SET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_LABEL_SET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1,"ModelUserNumber", &
                      & ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                  ENDIF
                  !Set pointers
                  IF(CELLML%MODELS_FIELD%MODELS_FIELD_AUTO_CREATED) THEN            
                    MODELS_FIELD=>CELLML%MODELS_FIELD%MODELS_FIELD
                  ELSE
                    CELLML%MODELS_FIELD%MODELS_FIELD=>MODELS_FIELD
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("CellML fields map source field domain is not associated.",ERR,ERROR,*999)         
                ENDIF
              ELSE
                CALL FLAG_ERROR("CellML fields map source geometric field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("CellML environment region is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("The CellML environment fields map has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment fields map is not associated. You must create the CellML field maps first.", &
          & ERR,ERROR,*999)
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
  SUBROUTINE CELLML_MODELS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we need to finish creation of the models field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
        IF(CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED) THEN
          CALL FLAG_ERROR("CellML models field has already been finished.",ERR,ERROR,*999)
        ELSE
          !Finish the models field creation
          IF(CELLML%MODELS_FIELD%MODELS_FIELD_AUTO_CREATED) &
            & CALL FIELD_CREATE_FINISH(CELLML%MODELS_FIELD%MODELS_FIELD,ERR,ERROR,*999)
          CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED=.TRUE.
          !Default the models field to the first model
          CALL FIELD_COMPONENT_VALUES_INITIALISE(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,1,1,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment models field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environement is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Finalise a CellML environment models field and deallocate all memory.
  SUBROUTINE CELLML_MODELS_FIELD_FINALISE(MODELS_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_MODELS_FIELD_TYPE), POINTER :: MODELS_FIELD !<A pointer to the CellML environment models field to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_MODELS_FIELD_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(MODELS_FIELD)) THEN
      DEALLOCATE(MODELS_FIELD)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_MODELS_FIELD_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_MODELS_FIELD_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_FINALISE

  !
  !=================================================================================================================================
  !

  !>Returns the models field for the given CellML environment.
  SUBROUTINE CELLML_MODELS_FIELD_GET(CELLML,MODELS_FIELD,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the models field.
    TYPE(FIELD_TYPE), POINTER :: MODELS_FIELD !<On return, a pointer to the models field for this CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_MODELS_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
        IF(CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED) THEN
          IF(ASSOCIATED(MODELS_FIELD)) THEN
            CALL FLAG_ERROR("Models field is already associated.",ERR,ERROR,*999)
          ELSE
            MODELS_FIELD=>CELLML%MODELS_FIELD%MODELS_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment models field has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment models field is not associated. Create the models field first.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Initialise a CellML environment models field.
  SUBROUTINE CELLML_MODELS_FIELD_INITIALISE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to initialise the models field for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_MODELS_FIELD_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
        CALL FLAG_ERROR("CellML environment models field is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CELLML%MODELS_FIELD,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML environment models field.",ERR,ERROR,*999)
        CELLML%MODELS_FIELD%CELLML=>CELLML
        CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED=.FALSE.
        CELLML%MODELS_FIELD%MODELS_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(CELLML%MODELS_FIELD%MODELS_FIELD)
        CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX=CELLML_MODELS_FIELD_NOT_CHECKED
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*998)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_MODELS_FIELD_INITIALISE")
    RETURN
999 CALL CELLML_MODELS_FIELD_FINALISE(CELLML%MODELS_FIELD,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_MODELS_FIELD_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_MODELS_FIELD_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_MODELS_FIELD_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Start the creation of the state field for the given CellML environment.
  SUBROUTINE CELLML_STATE_FIELD_CREATE_START(STATE_FIELD_USER_NUMBER,CELLML,STATE_FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: STATE_FIELD_USER_NUMBER !<The unique identifier for the state field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to create the state field.
    TYPE(FIELD_TYPE), POINTER :: STATE_FIELD  !<If associated on entry, a pointer to the user created state field which has the same user number as the specified state field user number. If not associated on entry, on exit, a pointer to the created state field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: component_idx
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,STATE_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML
    
    IF(ASSOCIATED(CELLML)) THEN
      CELLML_FIELD_MAPS=>CELLML%FIELD_MAPS
      IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
        IF(CELLML_FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
          IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
            CALL FLAG_ERROR("The CellML environment models field is already associated.",ERR,ERROR,*999)
          ELSE
            REGION=>CELLML%REGION
            IF(ASSOCIATED(REGION)) THEN
              IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)) THEN
                IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN)) THEN
                  IF(ASSOCIATED(STATE_FIELD)) THEN
                    !Check the field has been finished
                    IF(STATE_FIELD%FIELD_FINISHED) THEN
                      !Check the user numbers match
                      IF(STATE_FIELD_USER_NUMBER/=STATE_FIELD%USER_NUMBER) THEN
                        LOCAL_ERROR="The specified state field user number of "// &
                          & TRIM(NUMBER_TO_VSTRING(STATE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the user number of the specified state field of "// &
                          & TRIM(NUMBER_TO_VSTRING(STATE_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      STATE_FIELD_REGION=>STATE_FIELD%REGION
                      IF(ASSOCIATED(STATE_FIELD_REGION)) THEN
                        !Check the field is defined on the same region as the CellML region
                        IF(STATE_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                          LOCAL_ERROR="Invalid region setup. The specified state field has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(STATE_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & " and the CellML environment has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                        !Check the specified models field has the same geometric field as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,STATE_FIELD%GEOMETRIC_FIELD)) THEN
                          CALL FLAG_ERROR("The specified state field does not have the same geometric field as the "// &
                            & "geometric field for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                        !Check the specified models field has the same decomposition as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,STATE_FIELD%DECOMPOSITION)) THEN
                          CALL FLAG_ERROR("The specified state field does not have the same decomposition as the source "// &
                            & "domain decomposition for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("The specified state field region is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("The specified state field has not been finished.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    !Check the user number has not already been used for a field in this region.
                    NULLIFY(FIELD)
                    CALL FIELD_USER_NUMBER_FIND(STATE_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
                    IF(ASSOCIATED(FIELD)) THEN
                      LOCAL_ERROR="The specified state field user number of "// &
                        & TRIM(NUMBER_TO_VSTRING(STATE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                        & "has already been used to create a field on region number "// &
                        & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                  CALL CELLML_STATE_FIELD_INITIALISE(CELLML,ERR,ERROR,*999)
                  IF(ASSOCIATED(STATE_FIELD)) THEN
                    !Now check the supplied field.
                    CALL FIELD_DATA_TYPE_CHECK(STATE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_TYPE_CHECK(STATE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_CHECK(STATE_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_CHECK(STATE_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(STATE_FIELD,FIELD_U_VARIABLE_TYPE,CELLML%MAXIMUM_NUMBER_OF_STATE, &
                      & ERR,ERROR,*999)
                    DO component_idx=1,CELLML%MAXIMUM_NUMBER_OF_STATE
                      CALL FIELD_COMPONENT_MESH_COMPONENT_CHECK(STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ELSE
                    CELLML%STATE_FIELD%STATE_FIELD_AUTO_CREATED=.TRUE.
                    !Create the CellML environment models field
                    CALL FIELD_CREATE_START(STATE_FIELD_USER_NUMBER,REGION,CELLML%STATE_FIELD%STATE_FIELD,ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                      & ERR,ERROR,*999)
                    CALL FIELD_LABEL_SET(CELLML%STATE_FIELD%STATE_FIELD,"CellMLStateField",ERR,ERROR,*999)
                    CALL FIELD_TYPE_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,ERR,ERROR,*999)
                    CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,CELLML_FIELD_MAPS% &
                      & SOURCE_GEOMETRIC_FIELD,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_LABEL_SET(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE,"StateVariable", &
                      & ERR,ERROR,*999)
                    CALL FIELD_DOF_ORDER_TYPE_SET(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE,&
                      & CELLML%MAXIMUM_NUMBER_OF_STATE,ERR,ERROR,*999)
                    DO component_idx=1,CELLML%MAXIMUM_NUMBER_OF_STATE
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ENDIF
                  !Set pointers
                  IF(CELLML%STATE_FIELD%STATE_FIELD_AUTO_CREATED) THEN            
                    STATE_FIELD=>CELLML%STATE_FIELD%STATE_FIELD
                  ELSE
                    CELLML%STATE_FIELD%STATE_FIELD=>STATE_FIELD
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("CellML field maps source field domain is not associated.",ERR,ERROR,*999)         
                ENDIF
              ELSE
                CALL FLAG_ERROR("CellML field mapssource geometric field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("CellML environment region is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("The CellML environment fields map has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment fields map is not associated. You must create the CellML field maps first.", &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated",ERR,ERROR,*999)
    ENDIF

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
    INTEGER(INTG) :: model_idx,models_dof_idx,state_component_idx,CELLML_VARIABLE_TYPE
    INTEGER(INTG), POINTER :: MODELS_DATA(:)
    REAL(DP) :: INITIAL_VALUE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: MODELS_VARIABLE
    TYPE(CELLML_MODEL_TYPE), POINTER :: MODEL
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_STATE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
        IF(CELLML%STATE_FIELD%STATE_FIELD_FINISHED) THEN
          CALL FLAG_ERROR("CellML state field has already been finished.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
            IF(CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED) THEN
              CALL CELLML_MODELS_FIELD_CHECK(CELLML%MODELS_FIELD,ERR,ERROR,*999)
              !Finish the state field creation
              IF(CELLML%STATE_FIELD%STATE_FIELD_AUTO_CREATED) &
                & CALL FIELD_CREATE_FINISH(CELLML%STATE_FIELD%STATE_FIELD,ERR,ERROR,*999)              
              !Set the default field values to the initial CellML values.
              IF(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
                !Only one model so optimise
                MODEL=>CELLML%MODELS(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX)%PTR
                IF(ASSOCIATED(MODEL)) THEN
                  DO state_component_idx=1,MODEL%NUMBER_OF_STATE
                    CELLML_VARIABLE_TYPE=MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE(CELLML_STATE_FIELD,ERR,ERROR)
                    ERR = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(MODEL%PTR,CELLML_VARIABLE_TYPE,&
                      & state_component_idx,INITIAL_VALUE)
                    IF(ERR /= 0) THEN
                      !problem getting the initial value
                      LOCAL_ERROR="Failed to get an initial value for state variable with index "//&
                        & TRIM(NUMBER_TO_VSTRING(state_component_idx,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !WRITE(*,*) '(single model) Initial value for state variable: ',state_component_idx,'; type: ',&
                    !  & CELLML_VARIABLE_TYPE,'; value = ',INITIAL_VALUE
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(CELLML%STATE_FIELD%STATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,state_component_idx,INITIAL_VALUE,ERR,ERROR,*999)
                  ENDDO !state_component_idx
                ELSE
                  LOCAL_ERROR="The model is not associated for model index "// &
                    & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                !Multiple models so go through each dof.
                IF(ASSOCIATED(CELLML%FIELD_MAPS)) THEN
                  NULLIFY(MODELS_VARIABLE)
                  CALL FIELD_VARIABLE_GET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,MODELS_VARIABLE,ERR,ERROR,*999)
                  NULLIFY(MODELS_DATA)
                  CALL FIELD_PARAMETER_SET_DATA_GET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,MODELS_DATA,ERR,ERROR,*999)
                  DO models_dof_idx=1,MODELS_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    model_idx=MODELS_DATA(models_dof_idx)
                    IF(model_idx>0) THEN
                      MODEL=>CELLML%MODELS(model_idx)%PTR
                      IF(ASSOCIATED(MODEL)) THEN
                        DO state_component_idx=1,MODEL%NUMBER_OF_STATE
                          CELLML_VARIABLE_TYPE=MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE(CELLML_STATE_FIELD,ERR,ERROR)
                          IF(ERR/=0) GOTO 999
                          ERR = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(MODEL%PTR,CELLML_VARIABLE_TYPE,&
                            & state_component_idx,INITIAL_VALUE)
                          IF(ERR /= 0) THEN
                            !problem getting the initial value
                            LOCAL_ERROR="Failed to get an initial value for state variable with index "//&
                              & TRIM(NUMBER_TO_VSTRING(state_component_idx,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                          !WRITE(*,*) '(multiple models) Initial value for state variable: ',state_component_idx,'; type: ',&
                          !  & CELLML_VARIABLE_TYPE,'; value = ',INITIAL_VALUE
                          CALL CellML_FieldModelDofSet(MODELS_VARIABLE,models_dof_idx,CELLML%STATE_FIELD%STATE_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,state_component_idx,INITIAL_VALUE,ERR,ERROR,*999)
                        ENDDO !state_component_idx
                      ELSE
                        LOCAL_ERROR="The model is not associated for model index "// &
                          & TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ENDDO !dofIdx
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,MODELS_DATA,ERR,ERROR,*999)
                ELSE
                    CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
                  ENDIF
              ENDIF
              CELLML%STATE_FIELD%STATE_FIELD_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("CellML environment models field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("CellML environment models field is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment state field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Finalise a CellML environment state field and deallocate all memory.
  SUBROUTINE CELLML_STATE_FIELD_FINALISE(STATE_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_STATE_FIELD_TYPE), POINTER :: STATE_FIELD !<A pointer to the CellML environment state field to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_STATE_FIELD_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(STATE_FIELD)) THEN
      DEALLOCATE(STATE_FIELD)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_STATE_FIELD_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_STATE_FIELD_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_FINALISE

  !
  !=================================================================================================================================
  !

  !>Returns the state field for the given CellML environment.
  SUBROUTINE CELLML_STATE_FIELD_GET(CELLML,STATE_FIELD,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the state field.
    TYPE(FIELD_TYPE), POINTER :: STATE_FIELD !<On successful return will be set to the state field for this CellML environment
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_STATE_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
        IF(CELLML%STATE_FIELD%STATE_FIELD_FINISHED) THEN
          IF(ASSOCIATED(STATE_FIELD)) THEN
            CALL FLAG_ERROR("State field is already associated.",ERR,ERROR,*999)
          ELSE
            STATE_FIELD=>CELLML%STATE_FIELD%STATE_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment state field has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment state field is not associated. Create the state field first.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Initialise a CellML environment models field.
  SUBROUTINE CELLML_STATE_FIELD_INITIALISE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to initialise the models field for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_STATE_FIELD_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%STATE_FIELD)) THEN
        CALL FLAG_ERROR("CellML environment state field is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CELLML%STATE_FIELD,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML environment state field.",ERR,ERROR,*999)
        CELLML%STATE_FIELD%CELLML=>CELLML
        CELLML%STATE_FIELD%STATE_FIELD_FINISHED=.FALSE.
        CELLML%STATE_FIELD%STATE_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(CELLML%STATE_FIELD%STATE_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*998)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_STATE_FIELD_INITIALISE")
    RETURN
999 CALL CELLML_STATE_FIELD_FINALISE(CELLML%STATE_FIELD,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_STATE_FIELD_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_STATE_FIELD_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_STATE_FIELD_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Find the component ID in the given field for the variable defined by the given variable ID in the provided CellML environment.
  !! This generic routine will be used to map variable ID's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model ID also needed?
  !! - because the CellML fields should all be set up to allow direct use in the CellML code, the component number matches the index of the given variable in its associated array in the CellML generated code.
 SUBROUTINE CELLML_FIELD_COMPONENT_GET_C(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,COMPONENT_USER_NUMBER,ERR,ERROR,*)
   !Argument variables
   TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
   INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model to map from.
   INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CMISS_CELLML
   CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The ID of the model variable which needs to be located in the provided field.
   INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On return, the field component for the model variable defined by the given ID.
   INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
   TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
   !Local variables
   INTEGER(INTG) :: CELLML_VARIABLE_INDEX
   INTEGER(C_INT) :: ERROR_C
   TYPE(CELLML_MODEL_TYPE), POINTER :: CELLML_MODEL
   TYPE(VARYING_STRING) :: LOCAL_ERROR
   CHARACTER(LEN=1,KIND=C_CHAR) :: C_NAME(MAXSTRLEN)
   !INTEGER(INTG) :: C_NAME_L

   CALL ENTERS("CELLML_FIELD_COMPONENT_GET_C",ERR,ERROR,*999)

#ifdef USECELLML

   IF(ASSOCIATED(CELLML)) THEN
     CELLML_MODEL=>CELLML%MODELS(MODEL_INDEX)%PTR
     IF(ASSOCIATED(CELLML_MODEL)) THEN
       CALL CMISSF2CString(VARIABLE_ID,C_NAME)
       ERROR_C = CELLML_MODEL_DEFINITION_GET_VARIABLE_INDEX(CELLML_MODEL%PTR,C_NAME,CELLML_VARIABLE_INDEX)
       IF(ERROR_C /= 0) THEN
         LOCAL_ERROR="Failed to get the index for CellML variable: "// &
           & VARIABLE_ID//"; with the error code: "//TRIM(NUMBER_TO_VSTRING(ERROR_C,"*",ERR,ERROR))
         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
       ENDIF
       COMPONENT_USER_NUMBER=CELLML_VARIABLE_INDEX
     ENDIF
   ELSE
     CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
   ENDIF

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

  !>Find the component ID in the given field for the variable defined by the given variable ID in the provided CellML environment.
  !! This generic routine will be used to map variable ID's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model ID also needed?
  SUBROUTINE CELLML_FIELD_COMPONENT_GET_VS(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,COMPONENT_USER_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model to map from.
    INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CMISS_CELLML
    TYPE(VARYING_STRING), INTENT(IN) :: VARIABLE_ID !<The ID of the model variable which needs to be located in the provided field.
    INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On return, the field component for the model variable defined by the given ID.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_FIELD_COMPONENT_GET_VS",ERR,ERROR,*999)

#ifdef USECELLML

    CALL CELLML_FIELD_COMPONENT_GET(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,CHAR(VARIABLE_ID),COMPONENT_USER_NUMBER,ERR,ERROR,*999)

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
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_CREATE_START(INTERMEDIATE_FIELD_USER_NUMBER,CELLML,INTERMEDIATE_FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: INTERMEDIATE_FIELD_USER_NUMBER !<The unique identifier for the intermediate field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
    TYPE(FIELD_TYPE), POINTER :: INTERMEDIATE_FIELD !<If associated on entry, a pointer to the user created intermediate field which has the same user number as the specified intermediate field user number. If not associated on entry, on exit, a pointer to the created intermediate field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: component_idx
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,INTERMEDIATE_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CELLML_FIELD_MAPS=>CELLML%FIELD_MAPS
      IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
        IF(CELLML_FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
          IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
            CALL FLAG_ERROR("The CellML environment models field is already associated.",ERR,ERROR,*999)
          ELSE
            REGION=>CELLML%REGION
            IF(ASSOCIATED(REGION)) THEN
              IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)) THEN
                IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN)) THEN
                  IF(ASSOCIATED(INTERMEDIATE_FIELD)) THEN
                    !Check the field has been finished
                    IF(INTERMEDIATE_FIELD%FIELD_FINISHED) THEN
                      !Check the user numbers match
                      IF(INTERMEDIATE_FIELD_USER_NUMBER/=INTERMEDIATE_FIELD%USER_NUMBER) THEN
                        LOCAL_ERROR="The specified intermediate field user number of "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERMEDIATE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the user number of the specified intermediate field of "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERMEDIATE_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      INTERMEDIATE_FIELD_REGION=>INTERMEDIATE_FIELD%REGION
                      IF(ASSOCIATED(INTERMEDIATE_FIELD_REGION)) THEN
                        !Check the field is defined on the same region as the CellML region
                        IF(INTERMEDIATE_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                          LOCAL_ERROR="Invalid region setup. The specified intermediate field has been created on region"// &
                            & " number "//TRIM(NUMBER_TO_VSTRING(INTERMEDIATE_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & " and the specified CellML environment has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                        !Check the specified intermediate field has the same geometric field as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,INTERMEDIATE_FIELD%GEOMETRIC_FIELD)) THEN
                          CALL FLAG_ERROR("The specified intermediate field does not have the same geometric field as the "// &
                            & "geometric field for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                        !Check the specified intermediate field has the same decomposition as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION, &
                          & INTERMEDIATE_FIELD%DECOMPOSITION)) THEN
                          CALL FLAG_ERROR("The specified intermediate field does not have the same decomposition as the source "// &
                            & "domain decomposition for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("The specified intermediate field region is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("The specified intermediate field has not been finished.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    !Check the user number has not already been used for a field in this region.
                    NULLIFY(FIELD)
                    CALL FIELD_USER_NUMBER_FIND(INTERMEDIATE_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
                    IF(ASSOCIATED(FIELD)) THEN
                      LOCAL_ERROR="The specified intermediate field user number of "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERMEDIATE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                        & "has already been used to create a field on region number "// &
                        & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                  CALL CELLML_INTERMEDIATE_FIELD_INITIALISE(CELLML,ERR,ERROR,*999)
                  IF(ASSOCIATED(INTERMEDIATE_FIELD)) THEN
                    !Now check the supplied field.
                    CALL FIELD_DATA_TYPE_CHECK(INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_TYPE_CHECK(INTERMEDIATE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_CHECK(INTERMEDIATE_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_CHECK(INTERMEDIATE_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE,ERR,ERROR,*999)
                    DO component_idx=1,CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE
                      CALL FIELD_COMPONENT_MESH_COMPONENT_CHECK(INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                        & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                        & CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ELSE
                    CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_AUTO_CREATED=.TRUE.
                    !Create the CellML environment intermediate field
                    CALL FIELD_CREATE_START(INTERMEDIATE_FIELD_USER_NUMBER,REGION,CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                      & ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_LABEL_SET(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,"CellMLIntermediateField",ERR,ERROR,*999)
                    CALL FIELD_TYPE_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,ERR,ERROR,*999)
                    CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,CELLML_FIELD_MAPS% &
                      & SOURCE_GEOMETRIC_FIELD,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,[FIELD_U_VARIABLE_TYPE], &
                      & ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_LABEL_SET(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & "IntermediateVariable",ERR,ERROR,*999)
                    CALL FIELD_DOF_ORDER_TYPE_SET(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE,ERR,ERROR,*999)
                    DO component_idx=1,CELLML%MAXIMUM_NUMBER_OF_INTERMEDIATE
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ENDIF
                  !Set pointers
                  IF(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_AUTO_CREATED) THEN            
                    INTERMEDIATE_FIELD=>CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD
                  ELSE
                    CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD=>INTERMEDIATE_FIELD
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("CellML field maps source field domain is not associated.",ERR,ERROR,*999)         
                ENDIF
              ELSE
                CALL FLAG_ERROR("CellML field maps source geometric field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("CellML environment region is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("The CellML environment fields map has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment fields map is not associated. You must create the CellML field maps first.", &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated",ERR,ERROR,*999)
    ENDIF
    
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

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
        IF(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_FINISHED) THEN
          CALL FLAG_ERROR("CellML intermediate field has already been finished.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
            IF(CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED) THEN
              CALL CELLML_MODELS_FIELD_CHECK(CELLML%MODELS_FIELD,ERR,ERROR,*999)
              !Finish the intermediate field creation
              IF(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_AUTO_CREATED) &
                & CALL FIELD_CREATE_FINISH(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD,ERR,ERROR,*999)
              !As the intermediate field is strictly output do not initialise the values.
              CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("CellML environment models field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("CellML environment models field is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment intermediate field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Finalise a CellML environment models field and deallocate all memory.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_FINALISE(INTERMEDIATE_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_INTERMEDIATE_FIELD_TYPE), POINTER :: INTERMEDIATE_FIELD !<A pointer to the CellML environment intermediate field to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(INTERMEDIATE_FIELD)) THEN
      DEALLOCATE(INTERMEDIATE_FIELD)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_FINALISE

   !
  !=================================================================================================================================
  !

  !>Returns the intermediate field for the given CellML environment.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_GET(CELLML,INTERMEDIATE_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the intermediate field.
    TYPE(FIELD_TYPE), POINTER :: INTERMEDIATE_FIELD !<On return, a pointer to the intermediate field for this CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
        IF(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_FINISHED) THEN
          IF(ASSOCIATED(INTERMEDIATE_FIELD)) THEN
            CALL FLAG_ERROR("Intermediate field is already associated.",ERR,ERROR,*999)
          ELSE
            INTERMEDIATE_FIELD=>CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment intermediate field has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment intermediate field is not associated. Create the intermediate field first.", &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF
    
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

  !>Initialise a CellML environment intermediate field.
  SUBROUTINE CELLML_INTERMEDIATE_FIELD_INITIALISE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to initialise the intermediate field for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_INTERMEDIATE_FIELD_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%INTERMEDIATE_FIELD)) THEN
        CALL FLAG_ERROR("CellML environment intermediate field is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CELLML%INTERMEDIATE_FIELD,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML environment intermediate field.",ERR,ERROR,*999)
        CELLML%INTERMEDIATE_FIELD%CELLML=>CELLML
        CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_FINISHED=.FALSE.
        CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(CELLML%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*998)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_INTERMEDIATE_FIELD_INITIALISE")
    RETURN
999 CALL CELLML_INTERMEDIATE_FIELD_FINALISE(CELLML%INTERMEDIATE_FIELD,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_INTERMEDIATE_FIELD_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_INTERMEDIATE_FIELD_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_INTERMEDIATE_FIELD_INITIALISE

  !
  !=================================================================================================================================
  !

  !>Start the creation of the parameters field for the given CellML environment.
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_START(PARAMETERS_FIELD_USER_NUMBER,CELLML,PARAMETERS_FIELD,ERR,ERROR,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARAMETERS_FIELD_USER_NUMBER !<The unique identifier for the parameters field to be created for the given CellML environment object.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which we will be defining the parameters field.
    TYPE(FIELD_TYPE), POINTER :: PARAMETERS_FIELD  !<If associated on entry, a pointer to the user created parameters field which has the same user number as the specified parameters field user number. If not associated on entry, on exit, a pointer to the created parameters field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: component_idx
    TYPE(CELLML_FIELD_MAPS_TYPE), POINTER :: CELLML_FIELD_MAPS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,PARAMETERS_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_START",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      CELLML_FIELD_MAPS=>CELLML%FIELD_MAPS
      IF(ASSOCIATED(CELLML_FIELD_MAPS)) THEN
        IF(CELLML_FIELD_MAPS%CELLML_FIELD_MAPS_FINISHED) THEN
          IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
            CALL FLAG_ERROR("The CellML environment parameters field is already associated.",ERR,ERROR,*999)
          ELSE
            REGION=>CELLML%REGION
            IF(ASSOCIATED(REGION)) THEN
              IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD)) THEN
                IF(ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN)) THEN
                  IF(ASSOCIATED(PARAMETERS_FIELD)) THEN
                    !Check the field has been finished
                    IF(PARAMETERS_FIELD%FIELD_FINISHED) THEN
                      !Check the user numbers match
                      IF(PARAMETERS_FIELD_USER_NUMBER/=PARAMETERS_FIELD%USER_NUMBER) THEN
                        LOCAL_ERROR="The specified parameters field user number of "// &
                          & TRIM(NUMBER_TO_VSTRING(PARAMETERS_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not match the user number of the specified parameters field of "// &
                          & TRIM(NUMBER_TO_VSTRING(PARAMETERS_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      PARAMETERS_FIELD_REGION=>PARAMETERS_FIELD%REGION
                      IF(ASSOCIATED(PARAMETERS_FIELD_REGION)) THEN
                        !Check the field is defined on the same region as the CellML region
                        IF(PARAMETERS_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                          LOCAL_ERROR="Invalid region setup. The specified parameters field has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(PARAMETERS_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & " and the specified CellML environment has been created on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                        !Check the specified parameters field has the same geometric field as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,PARAMETERS_FIELD%GEOMETRIC_FIELD)) THEN
                          CALL FLAG_ERROR("The specified parameters field does not have the same geometric field as the "// &
                            & "geometric field for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                        !Check the specified parameters field has the same decomposition as the source field
                        IF(.NOT.ASSOCIATED(CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,PARAMETERS_FIELD%DECOMPOSITION)) THEN
                          CALL FLAG_ERROR("The specified parameters field does not have the same decomposition as the source "// &
                            & "domain decomposition for the specified CellML environment.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("The specified parameters field region is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("The specified parameters field has not been finished.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    !Check the user number has not already been used for a field in this region.
                    NULLIFY(FIELD)
                    CALL FIELD_USER_NUMBER_FIND(PARAMETERS_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
                    IF(ASSOCIATED(FIELD)) THEN
                      LOCAL_ERROR="The specified parameters field user number of "// &
                        & TRIM(NUMBER_TO_VSTRING(PARAMETERS_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                        & "has already been used to create a field on region number "// &
                        & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                  CALL CELLML_PARAMETERS_FIELD_INITIALISE(CELLML,ERR,ERROR,*999)
                  IF(ASSOCIATED(PARAMETERS_FIELD)) THEN
                    !Now check the supplied field.
                    CALL FIELD_DATA_TYPE_CHECK(PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,ERR,ERROR,*999)
                    CALL FIELD_TYPE_CHECK(PARAMETERS_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_CHECK(PARAMETERS_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_CHECK(PARAMETERS_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & CELLML%MAXIMUM_NUMBER_OF_PARAMETERS,ERR,ERROR,*999)
                    DO component_idx=1,CELLML%MAXIMUM_NUMBER_OF_PARAMETERS
                      CALL FIELD_COMPONENT_MESH_COMPONENT_CHECK(PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                        & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                        & CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ELSE
                    CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_AUTO_CREATED=.TRUE.
                    !Create the CellML environment parameters field
                    CALL FIELD_CREATE_START(PARAMETERS_FIELD_USER_NUMBER,REGION,CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                      & ERR,ERROR,*999)
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_DP_TYPE,ERR,ERROR,*999)
                    CALL FIELD_LABEL_SET(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,"CellMLParametersField",ERR,ERROR,*999)
                    CALL FIELD_TYPE_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                      & CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%DECOMPOSITION,ERR,ERROR,*999)
                    CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                      & CELLML_FIELD_MAPS%SOURCE_GEOMETRIC_FIELD,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,1,ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                      & ERR,ERROR,*999)
                    CALL FIELD_VARIABLE_LABEL_SET(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & "ParametersVariable",ERR,ERROR,*999)
                    CALL FIELD_DOF_ORDER_TYPE_SET(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & CELLML%MAXIMUM_NUMBER_OF_PARAMETERS,ERR,ERROR,*999)
                    DO component_idx=1,CELLML%MAXIMUM_NUMBER_OF_PARAMETERS
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_DOMAIN%MESH_COMPONENT_NUMBER, &
                        & ERR,ERROR,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                        & FIELD_U_VARIABLE_TYPE,component_idx,CELLML_FIELD_MAPS%SOURCE_FIELD_INTERPOLATION_TYPE,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ENDIF
                  !Set pointers
                  IF(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_AUTO_CREATED) THEN            
                    PARAMETERS_FIELD=>CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD
                  ELSE
                    CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD=>PARAMETERS_FIELD
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("CellML field maps source field domain is not associated.",ERR,ERROR,*999)         
                ENDIF
              ELSE
                CALL FLAG_ERROR("CellML field maps source geometric field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("CellML environment region is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("The CellML environment fields map has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment fields map is not associated. You must create the CellML field maps first.", &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated",ERR,ERROR,*999)
    ENDIF

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

  !>Finish the creation of the parameters field for the given CellML environment.
  SUBROUTINE CELLML_PARAMETERS_FIELD_CREATE_FINISH(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object for which to finalise the parameters field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables
    INTEGER(INTG) :: model_idx,models_dof_idx,parameter_component_idx,CELLML_VARIABLE_TYPE
    INTEGER(INTG), POINTER :: MODELS_DATA(:)
    REAL(DP) :: INITIAL_VALUE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: MODELS_VARIABLE
    TYPE(CELLML_MODEL_TYPE), POINTER :: MODEL
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CELLML_PARAMETERS_FIELD_CREATE_FINISH",ERR,ERROR,*999)

#ifdef USECELLML
    
    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
        IF(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_FINISHED) THEN
          CALL FLAG_ERROR("CellML parameters field has already been finished.",ERR,ERROR,*999)
        ELSE
         IF(ASSOCIATED(CELLML%MODELS_FIELD)) THEN
            IF(CELLML%MODELS_FIELD%MODELS_FIELD_FINISHED) THEN
              CALL CELLML_MODELS_FIELD_CHECK(CELLML%MODELS_FIELD,ERR,ERROR,*999)
              !Finish the parameters field creation
              IF(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_AUTO_CREATED) &
                & CALL FIELD_CREATE_FINISH(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,ERR,ERROR,*999)
              IF(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
                !Only one model so optimise
                MODEL=>CELLML%MODELS(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX)%PTR
                IF(ASSOCIATED(MODEL)) THEN
                  DO parameter_component_idx=1,MODEL%NUMBER_OF_PARAMETERS
                    CELLML_VARIABLE_TYPE=MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE(CELLML_PARAMETERS_FIELD,ERR,ERROR)
                    ERR = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(MODEL%PTR,CELLML_VARIABLE_TYPE,&
                      & parameter_component_idx,INITIAL_VALUE)
                    IF(ERR /= 0) THEN
                      !problem getting the initial value
                      LOCAL_ERROR="Failed to get an initial value for parameter variable with index "//&
                        & TRIM(NUMBER_TO_VSTRING(parameter_component_idx,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    !WRITE(*,*) '(single model) Initial value for parameter variable: ',parameter_component_idx,'; type: ',&
                    !  & CELLML_VARIABLE_TYPE,'; value = ',INITIAL_VALUE
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,parameter_component_idx,INITIAL_VALUE,ERR,ERROR,*999)
                  ENDDO !parameter_component_idx
                ELSE
                  LOCAL_ERROR="The model is not associated for model index "// &
                    & TRIM(NUMBER_TO_VSTRING(CELLML%MODELS_FIELD%ONLY_ONE_MODEL_INDEX,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                !Multiple models so go through each dof.
                IF(ASSOCIATED(CELLML%FIELD_MAPS)) THEN
                  NULLIFY(MODELS_VARIABLE)
                  CALL FIELD_VARIABLE_GET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE,MODELS_VARIABLE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_DATA_GET(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,MODELS_DATA,ERR,ERROR,*999)
                  DO models_dof_idx=1,MODELS_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    model_idx=MODELS_DATA(models_dof_idx)
                    IF(model_idx==0) THEN
                      ! Do nothing- empty model index specified
                    ELSE IF(model_idx > 0 .AND. model_idx <= CELLML%NUMBER_OF_MODELS) THEN
                      MODEL=>CELLML%MODELS(model_idx)%PTR
                      IF(ASSOCIATED(MODEL)) THEN
                        DO parameter_component_idx=1,MODEL%NUMBER_OF_PARAMETERS
                          CELLML_VARIABLE_TYPE=MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE(CELLML_PARAMETERS_FIELD,ERR,ERROR)
                          ERR = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(MODEL%PTR,CELLML_VARIABLE_TYPE,&
                            & parameter_component_idx,INITIAL_VALUE)
                          IF(ERR /= 0) THEN
                            !problem getting the initial value
                            LOCAL_ERROR="Failed to get an initial value for parameter variable with index "//&
                              & TRIM(NUMBER_TO_VSTRING(parameter_component_idx,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                          !WRITE(*,*) '(multiple models) Initial value for parameter variable: ',parameter_component_idx,'; type: ',&
                          !  & CELLML_VARIABLE_TYPE,'; value = ',INITIAL_VALUE
                          CALL CellML_FieldModelDofSet(MODELS_VARIABLE,models_dof_idx,CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parameter_component_idx,INITIAL_VALUE, &
                            & ERR,ERROR,*999)
                        ENDDO !parameter_component_idx
                      ELSE
                        LOCAL_ERROR="The model is not associated for model index "// &
                          & TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Invalid CellML model index: "// &
                        & TRIM(NUMBER_TO_VSTRING(model_idx,"*",ERR,ERROR))//". The specified index should be between 1 and "// &
                        & TRIM(NUMBER_TO_VSTRING(CELLML%NUMBER_OF_MODELS,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !models_dof_idx
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(CELLML%MODELS_FIELD%MODELS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,MODELS_DATA,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("CellML environment field maps is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("CellML environment models field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("CellML environment models field is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment parameters field is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environement is not associated.",ERR,ERROR,*999)
    ENDIF
    
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

  !>Finalise a CellML environment parameters field and deallocate all memory.
  SUBROUTINE CELLML_PARAMETERS_FIELD_FINALISE(PARAMETERS_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_PARAMETERS_FIELD_TYPE), POINTER :: PARAMETERS_FIELD !<A pointer to the CellML environment parameters field to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    
    CALL ENTERS("CELLML_PARAMETERS_FIELD_FINALISE",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(PARAMETERS_FIELD)) THEN
      DEALLOCATE(PARAMETERS_FIELD)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*999)

#endif

    CALL EXITS("CELLML_PARAMETERS_FIELD_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_PARAMETERS_FIELD_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_FINALISE

  !
  !=================================================================================================================================
  !

  !>Returns the parameters field for the given CellML environment.
  SUBROUTINE CELLML_PARAMETERS_FIELD_GET(CELLML,PARAMETERS_FIELD,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the parameters field.
    TYPE(FIELD_TYPE), POINTER :: PARAMETERS_FIELD !<On return, a pointer to the parameters field for this CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local variables

    CALL ENTERS("CELLML_PARAMETERS_FIELD_GET",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
        IF(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_FINISHED) THEN
          IF(ASSOCIATED(PARAMETERS_FIELD)) THEN
            CALL FLAG_ERROR("Parameters field is already associated.",ERR,ERROR,*999)
          ELSE
            PARAMETERS_FIELD=>CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("CellML environment parameters field has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("CellML environment parameters field is not associated. Create the parameters field first.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Initialise a CellML environment parameters field.
  SUBROUTINE CELLML_PARAMETERS_FIELD_INITIALISE(CELLML,ERR,ERROR,*)
    !Argument variables
    TYPE(CELLML_TYPE), POINTER :: CELLML !<A pointer to the CellML environment to initialise the models field for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !< The error string
    !Local variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CELLML_PARAMETERS_FIELD_INITIALISE",ERR,ERROR,*998)

#ifdef USECELLML

    IF(ASSOCIATED(CELLML)) THEN
      IF(ASSOCIATED(CELLML%PARAMETERS_FIELD)) THEN
        CALL FLAG_ERROR("CellML environment parameters field is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CELLML%PARAMETERS_FIELD,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate CellML environment parameters field.",ERR,ERROR,*999)
        CELLML%PARAMETERS_FIELD%CELLML=>CELLML
        CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_FINISHED=.FALSE.
        CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(CELLML%PARAMETERS_FIELD%PARAMETERS_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*998)
    ENDIF

#else

    CALL FLAG_ERROR("Must compile with USECELLML=true to use CellML functionality.",ERR,ERROR,*998)

#endif

    CALL EXITS("CELLML_PARAMETERS_FIELD_INITIALISE")
    RETURN
999 CALL CELLML_PARAMETERS_FIELD_FINALISE(CELLML%PARAMETERS_FIELD,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_PARAMETERS_FIELD_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_PARAMETERS_FIELD_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_PARAMETERS_FIELD_INITIALISE

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

    IF(ASSOCIATED(CELLML)) THEN
      IF(CELLML%CELLML_FINISHED) THEN 
        !Set the generated flag
        CELLML%CELLML_GENERATED=.TRUE.
      ELSE
        CALL FLAG_ERROR("CellML environment has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("CellML environment is not associated.",ERR,ERROR,*999)
    ENDIF

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

  !>Finds and returns in CELLML a pointer to the CellML environment identified by USER_NUMBER on a region. If no CellML environment with that USER_NUMBER exists CELLML is left nullified.
  SUBROUTINE CELLML_USER_NUMBER_FIND(USER_NUMBER,REGION,CELLML,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to find the CellML user number.
    TYPE(CELLML_TYPE), POINTER :: CELLML !<On return a pointer to the CellML environment with the given user number. If no CellML environment with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: cellml_idx
    TYPE(CELLML_ENVIRONMENTS_TYPE), POINTER :: CELLML_ENVIRONMENTS

    CALL ENTERS("CELLML_USER_NUMBER_FIND",ERR,ERROR,*999)

#ifdef USECELLML

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(CELLML)) THEN
        CALL FLAG_ERROR("CellML is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(CELLML)
        CELLML_ENVIRONMENTS=>REGION%CELLML_ENVIRONMENTS        
        IF(ASSOCIATED(CELLML_ENVIRONMENTS)) THEN
          cellml_idx=1
          DO WHILE(cellml_idx<=CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS.AND..NOT.ASSOCIATED(CELLML))
            IF(CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              CELLML=>CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR
            ELSE
              cellml_idx=cellml_idx+1
            ENDIF
          ENDDO
        ELSE
          CALL FLAG_ERROR("Region CellML environments is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
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

  !
  !=================================================================================================================================
  !

  !>Finalises the CellML environments and deallocates all memory.
  SUBROUTINE CELLML_ENVIRONMENTS_FINALISE(CELLML_ENVIRONMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(CELLML_ENVIRONMENTS_TYPE), POINTER :: CELLML_ENVIRONMENTS !<A pointer to the CellML environments to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: cellml_idx

    CALL ENTERS("CELLML_ENVIRONMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CELLML_ENVIRONMENTS)) THEN
      IF(ALLOCATED(CELLML_ENVIRONMENTS%ENVIRONMENTS)) THEN
        DO cellml_idx=1,SIZE(CELLML_ENVIRONMENTS%ENVIRONMENTS,1)
          CALL CELLML_FINALISE(CELLML_ENVIRONMENTS%ENVIRONMENTS(cellml_idx)%PTR,ERR,ERROR,*999)
        ENDDO !cellml_idx
        DEALLOCATE(CELLML_ENVIRONMENTS%ENVIRONMENTS)
      ENDIF
      DEALLOCATE(CELLML_ENVIRONMENTS)
    ENDIF

    CALL EXITS("CELLML_ENVIRONMENTS_FINALISE")
    RETURN
999 CALL ERRORS("CELLML_ENVIRONMENTS_FINALISE",ERR,ERROR)
    CALL EXITS("CELLML_ENVIRONMENTS_FINALISE")
    RETURN 1
  END SUBROUTINE CELLML_ENVIRONMENTS_FINALISE
  
  !
  !=================================================================================================================================
  !

  !>Initialises the CellML environments.
  SUBROUTINE CELLML_ENVIRONMENTS_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION  !<A pointer to the region to initialise the CellML environments for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("CELLML_ENVIRONMENTS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%CELLML_ENVIRONMENTS)) THEN
        CALL FLAG_ERROR("Region CellML environments is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(REGION%CELLML_ENVIRONMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region CellML environments.",ERR,ERROR,*999)
        REGION%CELLML_ENVIRONMENTS%REGION=>REGION
        REGION%CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("CELLML_ENVIRONMENTS_INITIALISE")
    RETURN
999 CALL CELLML_ENVIRONMENTS_FINALISE(REGION%CELLML_ENVIRONMENTS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CELLML_ENVIRONMENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("CELLML_ENVIRONMENTS_INITIALISE")
    RETURN 1
  END SUBROUTINE CELLML_ENVIRONMENTS_INITIALISE

  !
  !================================================================================================================================
  !

  !> Maps a CellML variable type to a CellML field type (\see CELLML_FieldTypes,CMISS_CELLML)
  FUNCTION MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG(CELLML_VARIABLE_TYPE,ERR,ERROR)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CELLML_VARIABLE_TYPE !<The CellML variable type to map.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ! CellML variable types in OpenCMISS(cellml)/CellMLModelDefinitionF.h:
    !   state=1; known=2; wanted=3.
    !   independent=4 - but untested and maybe not working yet.
    CALL ENTERS("MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG",ERR,ERROR,*999)
    
    SELECT CASE(CELLML_VARIABLE_TYPE)
    CASE(1)
      MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG=CELLML_STATE_FIELD
    CASE(2)
      MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG=CELLML_PARAMETERS_FIELD
    CASE(3)
      MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG=CELLML_INTERMEDIATE_FIELD
    CASE(4)
      LOCAL_ERROR="CellML variable type "//TRIM(NUMBER_TO_VSTRING(CELLML_VARIABLE_TYPE,"*",ERR,ERROR))//&
      & " (independent variable) support not yet implemented"
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="CellML variable type "//TRIM(NUMBER_TO_VSTRING(CELLML_VARIABLE_TYPE,"*",ERR,ERROR))//&
      & " is invalid or not implemented."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    
    CALL EXITS("MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG")
    RETURN
999 CALL ERRORS("MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG",ERR,ERROR)
    CALL EXITS("MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG")
    RETURN

  END FUNCTION MAP_CELLML_VARIABLE_TYPE_TO_FIELD_TYPE_INTG

  !
  !================================================================================================================================
  !

  !> Maps a CellML field type to a CellML variable type (\see CELLML_FieldTypes,CMISS_CELLML)
  FUNCTION MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG(CELLML_FIELD_TYPE,ERR,ERROR)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The CellML field type to map.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ! CellML variable types in OpenCMISS(cellml)/CellMLModelDefinitionF.h:
    !   state=1; known=2; wanted=3.
    !   independent=4 - but untested and maybe not working yet.
    CALL ENTERS("MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG",ERR,ERROR,*999)
    SELECT CASE(CELLML_FIELD_TYPE)
    CASE(CELLML_STATE_FIELD)
      MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG=1
    CASE(CELLML_PARAMETERS_FIELD)
      MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG=2
    CASE(CELLML_INTERMEDIATE_FIELD)
      MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG=3
    !CASE(4)
    !  LOCAL_ERROR="CellML variable type "//TRIM(NUMBER_TO_VSTRING(CELLML_VARIABLE_TYPE,"*",ERR,ERROR))//&
    !  & " (independent variable) support not yet implemented"
    !  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="CellML field type "//TRIM(NUMBER_TO_VSTRING(CELLML_FIELD_TYPE,"*",ERR,ERROR))//&
      & " is invalid or not implemented"
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    CALL EXITS("MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG")
    RETURN
999 CALL ERRORS("MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG",ERR,ERROR)
    CALL EXITS("MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG")
    RETURN

  END FUNCTION MAP_CELLML_FIELD_TYPE_TO_VARIABLE_TYPE_INTG

END MODULE CMISS_CELLML
