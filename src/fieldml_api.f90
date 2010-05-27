!> \file
!> $Id$
!> \author Caton Little
!> \brief This module contains the interface to the fieldml library
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

!> Interfaces the fieldml library
MODULE FIELDML_API

  USE ISO_C_BINDING

  IMPLICIT NONE

  INTEGER(C_INT), PARAMETER :: FML_INVALID_HANDLE = -1

  INTEGER(C_INT), PARAMETER :: FML_ERR_NO_ERROR = 0

  INTEGER(C_INT), PARAMETER :: FML_ERR_UNKNOWN_OBJECT = 1001

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_OBJECT = 1002

  INTEGER(C_INT), PARAMETER :: FML_ERR_INCOMPLETE_OBJECT = 1003

  INTEGER(C_INT), PARAMETER :: FML_ERR_MISCONFIGURED_OBJECT = 1004

  INTEGER(C_INT), PARAMETER :: FML_ERR_ACCESS_VIOLATION = 1005

  INTEGER(C_INT), PARAMETER :: FML_ERR_FILE_READ = 1006

  INTEGER(C_INT), PARAMETER :: FML_ERR_FILE_WRITE = 1007

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_1 = 1101

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_2 = 1102

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_3 = 1103

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_4 = 1104

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_5 = 1105

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_6 = 1106

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_7 = 1107

  INTEGER(C_INT), PARAMETER :: FML_ERR_INVALID_PARAMETER_8 = 1108

  INTEGER(C_INT), PARAMETER :: FML_ERR_UNSUPPORTED = 2000

  INTEGER(C_INT), PARAMETER :: BOUNDS_UNKNOWN = 0
  INTEGER(C_INT), PARAMETER :: BOUNDS_DISCRETE_CONTIGUOUS = 1
  INTEGER(C_INT), PARAMETER :: BOUNDS_DISCRETE_ARBITRARY = 2

  INTEGER(C_INT), PARAMETER :: TYPE_UNKNOWN = 0
  INTEGER(C_INT), PARAMETER :: TYPE_TEXT = 1
  INTEGER(C_INT), PARAMETER :: TYPE_LINES = 2

  INTEGER(C_INT), PARAMETER :: DESCRIPTION_UNKNOWN = 0
  INTEGER(C_INT), PARAMETER :: DESCRIPTION_SEMIDENSE = 1

  INTEGER(C_INT), PARAMETER :: LOCATION_UNKNOWN = 0
  INTEGER(C_INT), PARAMETER :: LOCATION_INLINE = 1
  INTEGER(C_INT), PARAMETER :: LOCATION_FILE = 2

  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN = 0
  INTEGER(C_INT), PARAMETER :: FHT_ENSEMBLE_DOMAIN = 1
  INTEGER(C_INT), PARAMETER :: FHT_CONTINUOUS_DOMAIN = 2
  INTEGER(C_INT), PARAMETER :: FHT_MESH_DOMAIN = 3
  INTEGER(C_INT), PARAMETER :: FHT_CONTINUOUS_IMPORT = 4
  INTEGER(C_INT), PARAMETER :: FHT_ENSEMBLE_PARAMETERS = 5
  INTEGER(C_INT), PARAMETER :: FHT_CONTINUOUS_PARAMETERS = 6
  INTEGER(C_INT), PARAMETER :: FHT_CONTINUOUS_PIECEWISE = 7
  INTEGER(C_INT), PARAMETER :: FHT_CONTINUOUS_AGGREGATE = 8
  INTEGER(C_INT), PARAMETER :: FHT_CONTINUOUS_VARIABLE = 9
  INTEGER(C_INT), PARAMETER :: FHT_ENSEMBLE_VARIABLE = 10
  INTEGER(C_INT), PARAMETER :: FHT_REMOTE_ENSEMBLE_DOMAIN = 11
  INTEGER(C_INT), PARAMETER :: FHT_REMOTE_CONTINUOUS_DOMAIN = 12
  INTEGER(C_INT), PARAMETER :: FHT_REMOTE_ENSEMBLE_EVALUATOR = 13
  INTEGER(C_INT), PARAMETER :: FHT_REMOTE_CONTINUOUS_EVALUATOR = 14
  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN_ENSEMBLE_DOMAIN = 15
  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN_CONTINUOUS_DOMAIN = 16
  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN_ENSEMBLE_EVALUATOR = 17
  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN_CONTINUOUS_EVALUATOR = 18
  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN_ENSEMBLE_SOURCE = 19
  INTEGER(C_INT), PARAMETER :: FHT_UNKNOWN_CONTINUOUS_SOURCE = 20

  INTERFACE
    FUNCTION Fieldml_CreateFromFile( filename ) &
      & BIND(C,NAME="Fieldml_CreateFromFile")
      USE TYPES
      USE ISO_C_BINDING
      CHARACTER(KIND=C_CHAR) :: filename(*)
      TYPE(C_PTR) :: Fieldml_CreateFromFile
    END FUNCTION Fieldml_CreateFromFile

    FUNCTION Fieldml_Create( ) &
      & BIND(C,NAME="Fieldml_Create")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR) :: Fieldml_Create
    END FUNCTION Fieldml_Create

    FUNCTION Fieldml_SetDebug( handle, debug ) &
      & BIND(C,NAME="Fieldml_SetDebug")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: debug
      INTEGER(C_INT) :: Fieldml_SetDebug
    END FUNCTION Fieldml_SetDebug

    FUNCTION Fieldml_GetLastError( handle ) &
      & BIND(C,NAME="Fieldml_GetLastError")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT) :: Fieldml_GetLastError
    END FUNCTION Fieldml_GetLastError

    FUNCTION Fieldml_WriteFile( handle, filename ) &
      & BIND(C,NAME="Fieldml_WriteFile")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: filename(*)
      INTEGER(C_INT) :: Fieldml_WriteFile
    END FUNCTION Fieldml_WriteFile

    FUNCTION Fieldml_Destroy( handle ) &
      & BIND(C,NAME="Fieldml_Destroy")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT) :: Fieldml_Destroy
    END FUNCTION Fieldml_Destroy

    FUNCTION Fieldml_GetErrorCount( handle ) &
      & BIND(C,NAME="Fieldml_GetErrorCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT) :: Fieldml_GetErrorCount
    END FUNCTION Fieldml_GetErrorCount

    FUNCTION Fieldml_CopyError( handle, errorIndex, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyError")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: errorIndex
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyError
    END FUNCTION Fieldml_CopyError

    FUNCTION Fieldml_GetObjectCount( handle, type ) &
      & BIND(C,NAME="Fieldml_GetObjectCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: type
      INTEGER(C_INT) :: Fieldml_GetObjectCount
    END FUNCTION Fieldml_GetObjectCount

    FUNCTION Fieldml_GetObject( handle, type, objectIndex ) &
      & BIND(C,NAME="Fieldml_GetObject")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: type
      INTEGER(C_INT), VALUE :: objectIndex
      INTEGER(C_INT) :: Fieldml_GetObject
    END FUNCTION Fieldml_GetObject

    FUNCTION Fieldml_GetObjectType( handle, object ) &
      & BIND(C,NAME="Fieldml_GetObjectType")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: object
      INTEGER(C_INT) :: Fieldml_GetObjectType
    END FUNCTION Fieldml_GetObjectType

    FUNCTION Fieldml_GetNamedObject( handle, name ) &
      & BIND(C,NAME="Fieldml_GetNamedObject")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT) :: Fieldml_GetNamedObject
    END FUNCTION Fieldml_GetNamedObject

    FUNCTION Fieldml_GetMarkupCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetMarkupCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetMarkupCount
    END FUNCTION Fieldml_GetMarkupCount

    FUNCTION Fieldml_ValidateObject( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_ValidateObject")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_ValidateObject
    END FUNCTION Fieldml_ValidateObject

    FUNCTION Fieldml_CopyMarkupAttribute( handle, objectHandle, markupIndex, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyMarkupAttribute")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: markupIndex
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyMarkupAttribute
    END FUNCTION Fieldml_CopyMarkupAttribute

    FUNCTION Fieldml_CopyMarkupValue( handle, objectHandle, markupIndex, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyMarkupValue")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: markupIndex
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyMarkupValue
    END FUNCTION Fieldml_CopyMarkupValue

    FUNCTION Fieldml_CopyMarkupAttributeValue( handle, objectHandle, attribute, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyMarkupAttributeValue")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      CHARACTER(KIND=C_CHAR) :: attribute(*)
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyMarkupAttributeValue
    END FUNCTION Fieldml_CopyMarkupAttributeValue

    FUNCTION Fieldml_SetMarkup( handle, objectHandle, attribute, value ) &
      & BIND(C,NAME="Fieldml_SetMarkup")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      CHARACTER(KIND=C_CHAR) :: attribute(*)
      CHARACTER(KIND=C_CHAR) :: value(*)
      INTEGER(C_INT) :: Fieldml_SetMarkup
    END FUNCTION Fieldml_SetMarkup

    FUNCTION Fieldml_GetDomainComponentEnsemble( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetDomainComponentEnsemble")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetDomainComponentEnsemble
    END FUNCTION Fieldml_GetDomainComponentEnsemble

    FUNCTION Fieldml_CreateEnsembleDomain( handle, name, componentHandle ) &
      & BIND(C,NAME="Fieldml_CreateEnsembleDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: componentHandle
      INTEGER(C_INT) :: Fieldml_CreateEnsembleDomain
    END FUNCTION Fieldml_CreateEnsembleDomain

    FUNCTION Fieldml_CreateContinuousDomain( handle, name, componentHandle ) &
      & BIND(C,NAME="Fieldml_CreateContinuousDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: componentHandle
      INTEGER(C_INT) :: Fieldml_CreateContinuousDomain
    END FUNCTION Fieldml_CreateContinuousDomain

    FUNCTION Fieldml_CreateMeshDomain( handle, name, xiEnsemble ) &
      & BIND(C,NAME="Fieldml_CreateMeshDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: xiEnsemble
      INTEGER(C_INT) :: Fieldml_CreateMeshDomain
    END FUNCTION Fieldml_CreateMeshDomain

    FUNCTION Fieldml_GetMeshXiDomain( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetMeshXiDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetMeshXiDomain
    END FUNCTION Fieldml_GetMeshXiDomain

    FUNCTION Fieldml_GetMeshElementDomain( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetMeshElementDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetMeshElementDomain
    END FUNCTION Fieldml_GetMeshElementDomain

    FUNCTION Fieldml_CopyMeshElementShape( handle, objectHandle, elementNumber, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyMeshElementShape")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: elementNumber
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyMeshElementShape
    END FUNCTION Fieldml_CopyMeshElementShape

    FUNCTION Fieldml_SetMeshDefaultShape( handle, mesh, shape ) &
      & BIND(C,NAME="Fieldml_SetMeshDefaultShape")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: mesh
      CHARACTER(KIND=C_CHAR) :: shape(*)
      INTEGER(C_INT) :: Fieldml_SetMeshDefaultShape
    END FUNCTION Fieldml_SetMeshDefaultShape

    FUNCTION Fieldml_SetMeshElementShape( handle, mesh, elementNumber, shape ) &
      & BIND(C,NAME="Fieldml_SetMeshElementShape")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: mesh
      INTEGER(C_INT), VALUE :: elementNumber
      CHARACTER(KIND=C_CHAR) :: shape(*)
      INTEGER(C_INT) :: Fieldml_SetMeshElementShape
    END FUNCTION Fieldml_SetMeshElementShape

    FUNCTION Fieldml_GetMeshConnectivityCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetMeshConnectivityCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetMeshConnectivityCount
    END FUNCTION Fieldml_GetMeshConnectivityCount

    FUNCTION Fieldml_GetMeshConnectivityDomain( handle, objectHandle, connectivityIndex ) &
      & BIND(C,NAME="Fieldml_GetMeshConnectivityDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: connectivityIndex
      INTEGER(C_INT) :: Fieldml_GetMeshConnectivityDomain
    END FUNCTION Fieldml_GetMeshConnectivityDomain

    FUNCTION Fieldml_GetMeshConnectivitySource( handle, objectHandle, connectivityIndex ) &
      & BIND(C,NAME="Fieldml_GetMeshConnectivitySource")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: connectivityIndex
      INTEGER(C_INT) :: Fieldml_GetMeshConnectivitySource
    END FUNCTION Fieldml_GetMeshConnectivitySource

    FUNCTION Fieldml_SetMeshConnectivity( handle, mesh, evaluator, pointDomain ) &
      & BIND(C,NAME="Fieldml_SetMeshConnectivity")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: mesh
      INTEGER(C_INT), VALUE :: evaluator
      INTEGER(C_INT), VALUE :: pointDomain
      INTEGER(C_INT) :: Fieldml_SetMeshConnectivity
    END FUNCTION Fieldml_SetMeshConnectivity

    FUNCTION Fieldml_GetDomainBoundsType( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetDomainBoundsType")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetDomainBoundsType
    END FUNCTION Fieldml_GetDomainBoundsType

    FUNCTION Fieldml_GetEnsembleDomainElementCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetEnsembleDomainElementCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetEnsembleDomainElementCount
    END FUNCTION Fieldml_GetEnsembleDomainElementCount

    FUNCTION Fieldml_GetEnsembleDomainElementNames( handle, objectHandle, array, arrayLength ) &
      & BIND(C,NAME="Fieldml_GetEnsembleDomainElementNames")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      TYPE(C_PTR), VALUE :: array
      INTEGER(C_INT), VALUE :: arrayLength
      INTEGER(C_INT) :: Fieldml_GetEnsembleDomainElementNames
    END FUNCTION Fieldml_GetEnsembleDomainElementNames

    FUNCTION Fieldml_GetContiguousBoundsCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetContiguousBoundsCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetContiguousBoundsCount
    END FUNCTION Fieldml_GetContiguousBoundsCount

    FUNCTION Fieldml_SetContiguousBoundsCount( handle, objectHandle, count ) &
      & BIND(C,NAME="Fieldml_SetContiguousBoundsCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: count
      INTEGER(C_INT) :: Fieldml_SetContiguousBoundsCount
    END FUNCTION Fieldml_SetContiguousBoundsCount

    FUNCTION Fieldml_CopyObjectName( handle, objectHandle, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyObjectName")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyObjectName
    END FUNCTION Fieldml_CopyObjectName

    FUNCTION Fieldml_GetValueDomain( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetValueDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetValueDomain
    END FUNCTION Fieldml_GetValueDomain

    FUNCTION Fieldml_CreateEnsembleVariable( handle, name, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateEnsembleVariable")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateEnsembleVariable
    END FUNCTION Fieldml_CreateEnsembleVariable

    FUNCTION Fieldml_CreateContinuousVariable( handle, name, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateContinuousVariable")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateContinuousVariable
    END FUNCTION Fieldml_CreateContinuousVariable

    FUNCTION Fieldml_CreateEnsembleParameters( handle, name, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateEnsembleParameters")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateEnsembleParameters
    END FUNCTION Fieldml_CreateEnsembleParameters

    FUNCTION Fieldml_CreateContinuousParameters( handle, name, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateContinuousParameters")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateContinuousParameters
    END FUNCTION Fieldml_CreateContinuousParameters

    FUNCTION Fieldml_SetParameterDataDescription( handle, objectHandle, description ) &
      & BIND(C,NAME="Fieldml_SetParameterDataDescription")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: description
      INTEGER(C_INT) :: Fieldml_SetParameterDataDescription
    END FUNCTION Fieldml_SetParameterDataDescription

    FUNCTION Fieldml_GetParameterDataDescription( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetParameterDataDescription")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetParameterDataDescription
    END FUNCTION Fieldml_GetParameterDataDescription

    FUNCTION Fieldml_GetParameterDataLocation( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetParameterDataLocation")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetParameterDataLocation
    END FUNCTION Fieldml_GetParameterDataLocation

    FUNCTION Fieldml_SetParameterDataLocation( handle, objectHandle, location ) &
      & BIND(C,NAME="Fieldml_SetParameterDataLocation")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: location
      INTEGER(C_INT) :: Fieldml_SetParameterDataLocation
    END FUNCTION Fieldml_SetParameterDataLocation

    FUNCTION Fieldml_AddInlineParameterData( handle, objectHandle, data, length ) &
      & BIND(C,NAME="Fieldml_AddInlineParameterData")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      CHARACTER(KIND=C_CHAR) :: data(*)
      INTEGER(C_INT), VALUE :: length
      INTEGER(C_INT) :: Fieldml_AddInlineParameterData
    END FUNCTION Fieldml_AddInlineParameterData

    FUNCTION Fieldml_SetParameterFileData( handle, objectHandle, filename, type, offset ) &
      & BIND(C,NAME="Fieldml_SetParameterFileData")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      CHARACTER(KIND=C_CHAR) :: filename(*)
      INTEGER(C_INT), VALUE :: type
      INTEGER(C_INT), VALUE :: offset
      INTEGER(C_INT) :: Fieldml_SetParameterFileData
    END FUNCTION Fieldml_SetParameterFileData

    FUNCTION Fieldml_CopyParameterDataFilename( handle, objectHandle, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopyParameterDataFilename")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      CHARACTER(KIND=C_CHAR) :: buffer(*)
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopyParameterDataFilename
    END FUNCTION Fieldml_CopyParameterDataFilename

    FUNCTION Fieldml_GetParameterDataOffset( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetParameterDataOffset")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetParameterDataOffset
    END FUNCTION Fieldml_GetParameterDataOffset

    FUNCTION Fieldml_GetParameterDataFileType( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetParameterDataFileType")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetParameterDataFileType
    END FUNCTION Fieldml_GetParameterDataFileType

    FUNCTION Fieldml_AddSemidenseIndex( handle, objectHandle, indexHandle, isSparse ) &
      & BIND(C,NAME="Fieldml_AddSemidenseIndex")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: indexHandle
      INTEGER(C_INT), VALUE :: isSparse
      INTEGER(C_INT) :: Fieldml_AddSemidenseIndex
    END FUNCTION Fieldml_AddSemidenseIndex

    FUNCTION Fieldml_GetSemidenseIndexCount( handle, objectHandle, isSparse ) &
      & BIND(C,NAME="Fieldml_GetSemidenseIndexCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: isSparse
      INTEGER(C_INT) :: Fieldml_GetSemidenseIndexCount
    END FUNCTION Fieldml_GetSemidenseIndexCount

    FUNCTION Fieldml_GetSemidenseIndex( handle, objectHandle, indexIndex, isSparse ) &
      & BIND(C,NAME="Fieldml_GetSemidenseIndex")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: indexIndex
      INTEGER(C_INT), VALUE :: isSparse
      INTEGER(C_INT) :: Fieldml_GetSemidenseIndex
    END FUNCTION Fieldml_GetSemidenseIndex

    FUNCTION Fieldml_SetSwizzle( handle, objectHandle, buffer, count ) &
      & BIND(C,NAME="Fieldml_SetSwizzle")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      TYPE(C_PTR), VALUE :: buffer
      INTEGER(C_INT), VALUE :: count
      INTEGER(C_INT) :: Fieldml_SetSwizzle
    END FUNCTION Fieldml_SetSwizzle

    FUNCTION Fieldml_GetSwizzleCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetSwizzleCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetSwizzleCount
    END FUNCTION Fieldml_GetSwizzleCount

    FUNCTION Fieldml_CopySwizzleData( handle, objectHandle, buffer, bufferLength ) &
      & BIND(C,NAME="Fieldml_CopySwizzleData")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      TYPE(C_PTR), VALUE :: buffer
      INTEGER(C_INT), VALUE :: bufferLength
      INTEGER(C_INT) :: Fieldml_CopySwizzleData
    END FUNCTION Fieldml_CopySwizzleData

    FUNCTION Fieldml_CreateContinuousPiecewise( handle, name, indexHandle, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateContinuousPiecewise")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: indexHandle
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateContinuousPiecewise
    END FUNCTION Fieldml_CreateContinuousPiecewise

    FUNCTION Fieldml_CreateContinuousAggregate( handle, name, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateContinuousAggregate")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateContinuousAggregate
    END FUNCTION Fieldml_CreateContinuousAggregate

    FUNCTION Fieldml_SetDefaultEvaluator( handle, objectHandle, evaluator ) &
      & BIND(C,NAME="Fieldml_SetDefaultEvaluator")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: evaluator
      INTEGER(C_INT) :: Fieldml_SetDefaultEvaluator
    END FUNCTION Fieldml_SetDefaultEvaluator

    FUNCTION Fieldml_SetEvaluator( handle, objectHandle, element, evaluator ) &
      & BIND(C,NAME="Fieldml_SetEvaluator")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: element
      INTEGER(C_INT), VALUE :: evaluator
      INTEGER(C_INT) :: Fieldml_SetEvaluator
    END FUNCTION Fieldml_SetEvaluator

    FUNCTION Fieldml_GetEvaluatorCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetEvaluatorCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetEvaluatorCount
    END FUNCTION Fieldml_GetEvaluatorCount

    FUNCTION Fieldml_GetEvaluatorElement( handle, objectHandle, evaluatorIndex ) &
      & BIND(C,NAME="Fieldml_GetEvaluatorElement")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: evaluatorIndex
      INTEGER(C_INT) :: Fieldml_GetEvaluatorElement
    END FUNCTION Fieldml_GetEvaluatorElement

    FUNCTION Fieldml_GetEvaluator( handle, objectHandle, evaluatorIndex ) &
      & BIND(C,NAME="Fieldml_GetEvaluator")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: evaluatorIndex
      INTEGER(C_INT) :: Fieldml_GetEvaluator
    END FUNCTION Fieldml_GetEvaluator

    FUNCTION Fieldml_GetElementEvaluator( handle, objectHandle, elementNumber ) &
      & BIND(C,NAME="Fieldml_GetElementEvaluator")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: elementNumber
      INTEGER(C_INT) :: Fieldml_GetElementEvaluator
    END FUNCTION Fieldml_GetElementEvaluator

    FUNCTION Fieldml_CreateContinuousImport( handle, name, remoteEvaluator, valueDomain ) &
      & BIND(C,NAME="Fieldml_CreateContinuousImport")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(KIND=C_CHAR) :: name(*)
      INTEGER(C_INT), VALUE :: remoteEvaluator
      INTEGER(C_INT), VALUE :: valueDomain
      INTEGER(C_INT) :: Fieldml_CreateContinuousImport
    END FUNCTION Fieldml_CreateContinuousImport

    FUNCTION Fieldml_GetImportRemoteEvaluator( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetImportRemoteEvaluator")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetImportRemoteEvaluator
    END FUNCTION Fieldml_GetImportRemoteEvaluator

    FUNCTION Fieldml_GetAliasCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetAliasCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetAliasCount
    END FUNCTION Fieldml_GetAliasCount

    FUNCTION Fieldml_GetAliasLocal( handle, objectHandle, aliasIndex ) &
      & BIND(C,NAME="Fieldml_GetAliasLocal")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: aliasIndex
      INTEGER(C_INT) :: Fieldml_GetAliasLocal
    END FUNCTION Fieldml_GetAliasLocal

    FUNCTION Fieldml_GetAliasRemote( handle, objectHandle, aliasIndex ) &
      & BIND(C,NAME="Fieldml_GetAliasRemote")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: aliasIndex
      INTEGER(C_INT) :: Fieldml_GetAliasRemote
    END FUNCTION Fieldml_GetAliasRemote

    FUNCTION Fieldml_SetAlias( handle, objectHandle, remoteDomain, localSource ) &
      & BIND(C,NAME="Fieldml_SetAlias")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: remoteDomain
      INTEGER(C_INT), VALUE :: localSource
      INTEGER(C_INT) :: Fieldml_SetAlias
    END FUNCTION Fieldml_SetAlias

    FUNCTION Fieldml_GetIndexCount( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_GetIndexCount")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT) :: Fieldml_GetIndexCount
    END FUNCTION Fieldml_GetIndexCount

    FUNCTION Fieldml_GetIndexDomain( handle, objectHandle, indexIndex ) &
      & BIND(C,NAME="Fieldml_GetIndexDomain")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: indexIndex
      INTEGER(C_INT) :: Fieldml_GetIndexDomain
    END FUNCTION Fieldml_GetIndexDomain

    FUNCTION Fieldml_OpenReader( handle, objectHandle ) &
      & BIND(C,NAME="Fieldml_OpenReader")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      TYPE(C_PTR) :: Fieldml_OpenReader
    END FUNCTION Fieldml_OpenReader

    FUNCTION Fieldml_ReadIntSlice( handle, reader, indexBuffer, valueBuffer ) &
      & BIND(C,NAME="Fieldml_ReadIntSlice")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      TYPE(C_PTR), VALUE :: reader
      TYPE(C_PTR), VALUE :: indexBuffer
      TYPE(C_PTR), VALUE :: valueBuffer
      INTEGER(C_INT) :: Fieldml_ReadIntSlice
    END FUNCTION Fieldml_ReadIntSlice

    FUNCTION Fieldml_ReadDoubleSlice( handle, reader, indexBuffer, valueBuffer ) &
      & BIND(C,NAME="Fieldml_ReadDoubleSlice")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      TYPE(C_PTR), VALUE :: reader
      TYPE(C_PTR), VALUE :: indexBuffer
      TYPE(C_PTR), VALUE :: valueBuffer
      INTEGER(C_INT) :: Fieldml_ReadDoubleSlice
    END FUNCTION Fieldml_ReadDoubleSlice

    FUNCTION Fieldml_CloseReader( handle, reader ) &
      & BIND(C,NAME="Fieldml_CloseReader")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      TYPE(C_PTR), VALUE :: reader
      INTEGER(C_INT) :: Fieldml_CloseReader
    END FUNCTION Fieldml_CloseReader

    FUNCTION Fieldml_OpenWriter( handle, objectHandle, append ) &
      & BIND(C,NAME="Fieldml_OpenWriter")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT), VALUE :: objectHandle
      INTEGER(C_INT), VALUE :: append
      TYPE(C_PTR) :: Fieldml_OpenWriter
    END FUNCTION Fieldml_OpenWriter

    FUNCTION Fieldml_WriteIntSlice( handle, writer, indexBuffer, valueBuffer ) &
      & BIND(C,NAME="Fieldml_WriteIntSlice")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      TYPE(C_PTR), VALUE :: writer
      TYPE(C_PTR), VALUE :: indexBuffer
      TYPE(C_PTR), VALUE :: valueBuffer
      INTEGER(C_INT) :: Fieldml_WriteIntSlice
    END FUNCTION Fieldml_WriteIntSlice

    FUNCTION Fieldml_WriteDoubleSlice( handle, writer, indexBuffer, valueBuffer ) &
      & BIND(C,NAME="Fieldml_WriteDoubleSlice")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      TYPE(C_PTR), VALUE :: writer
      TYPE(C_PTR), VALUE :: indexBuffer
      TYPE(C_PTR), VALUE :: valueBuffer
      INTEGER(C_INT) :: Fieldml_WriteDoubleSlice
    END FUNCTION Fieldml_WriteDoubleSlice

    FUNCTION Fieldml_CloseWriter( handle, writer ) &
      & BIND(C,NAME="Fieldml_CloseWriter")
      USE TYPES
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      TYPE(C_PTR), VALUE :: writer
      INTEGER(C_INT) :: Fieldml_CloseWriter
    END FUNCTION Fieldml_CloseWriter

  END INTERFACE

  PUBLIC Fieldml_CreateFromFile, Fieldml_Create, Fieldml_SetDebug, Fieldml_GetLastError, Fieldml_WriteFile, &
    & Fieldml_Destroy, Fieldml_GetErrorCount, Fieldml_CopyError, Fieldml_GetObjectCount, Fieldml_GetObject, &
    & Fieldml_GetObjectType, Fieldml_GetNamedObject, Fieldml_GetMarkupCount, Fieldml_ValidateObject, &
    & Fieldml_CopyMarkupAttribute, Fieldml_CopyMarkupValue, Fieldml_CopyMarkupAttributeValue, Fieldml_SetMarkup, &
    & Fieldml_GetDomainComponentEnsemble, Fieldml_CreateEnsembleDomain, Fieldml_CreateContinuousDomain, &
    & Fieldml_CreateMeshDomain, Fieldml_GetMeshXiDomain, Fieldml_GetMeshElementDomain, Fieldml_CopyMeshElementShape, &
    & Fieldml_SetMeshDefaultShape, Fieldml_SetMeshElementShape, Fieldml_GetMeshConnectivityCount, &
    & Fieldml_GetMeshConnectivityDomain, Fieldml_GetMeshConnectivitySource, Fieldml_SetMeshConnectivity, &
    & Fieldml_GetDomainBoundsType, Fieldml_GetEnsembleDomainElementCount, Fieldml_GetEnsembleDomainElementNames, &
    & Fieldml_GetContiguousBoundsCount, Fieldml_SetContiguousBoundsCount, Fieldml_CopyObjectName, Fieldml_GetValueDomain, &
    & Fieldml_CreateEnsembleVariable, Fieldml_CreateContinuousVariable, Fieldml_CreateEnsembleParameters, &
    & Fieldml_CreateContinuousParameters, Fieldml_SetParameterDataDescription, Fieldml_GetParameterDataDescription, &
    & Fieldml_GetParameterDataLocation, Fieldml_SetParameterDataLocation, Fieldml_AddInlineParameterData, &
    & Fieldml_SetParameterFileData, Fieldml_CopyParameterDataFilename, Fieldml_GetParameterDataOffset, &
    & Fieldml_GetParameterDataFileType, Fieldml_AddSemidenseIndex, Fieldml_GetSemidenseIndexCount, &
    & Fieldml_GetSemidenseIndex, Fieldml_SetSwizzle, Fieldml_GetSwizzleCount, Fieldml_CopySwizzleData, &
    & Fieldml_CreateContinuousPiecewise, Fieldml_CreateContinuousAggregate, Fieldml_SetDefaultEvaluator, &
    & Fieldml_SetEvaluator, Fieldml_GetEvaluatorCount, Fieldml_GetEvaluatorElement, Fieldml_GetEvaluator, &
    & Fieldml_GetElementEvaluator, Fieldml_CreateContinuousImport, Fieldml_GetImportRemoteEvaluator, Fieldml_GetAliasCount, &
    & Fieldml_GetAliasLocal, Fieldml_GetAliasRemote, Fieldml_SetAlias, Fieldml_GetIndexCount, Fieldml_GetIndexDomain, &
    & Fieldml_OpenReader, Fieldml_ReadIntSlice, Fieldml_ReadDoubleSlice, Fieldml_CloseReader, Fieldml_OpenWriter, &
    & Fieldml_WriteIntSlice, Fieldml_WriteDoubleSlice, Fieldml_CloseWriter

  PUBLIC BOUNDS_UNKNOWN, BOUNDS_DISCRETE_CONTIGUOUS, BOUNDS_DISCRETE_ARBITRARY

  PUBLIC TYPE_UNKNOWN, TYPE_TEXT, TYPE_LINES

  PUBLIC DESCRIPTION_UNKNOWN, DESCRIPTION_SEMIDENSE

  PUBLIC LOCATION_UNKNOWN, LOCATION_INLINE, LOCATION_FILE

  PUBLIC FHT_UNKNOWN, FHT_ENSEMBLE_DOMAIN, FHT_CONTINUOUS_DOMAIN, FHT_MESH_DOMAIN, FHT_CONTINUOUS_IMPORT, &
    & FHT_ENSEMBLE_PARAMETERS, FHT_CONTINUOUS_PARAMETERS, FHT_CONTINUOUS_PIECEWISE, FHT_CONTINUOUS_AGGREGATE, &
    & FHT_CONTINUOUS_VARIABLE, FHT_ENSEMBLE_VARIABLE, FHT_REMOTE_ENSEMBLE_DOMAIN, FHT_REMOTE_CONTINUOUS_DOMAIN, &
    & FHT_REMOTE_ENSEMBLE_EVALUATOR, FHT_REMOTE_CONTINUOUS_EVALUATOR, FHT_UNKNOWN_ENSEMBLE_DOMAIN, &
    & FHT_UNKNOWN_CONTINUOUS_DOMAIN, FHT_UNKNOWN_ENSEMBLE_EVALUATOR, FHT_UNKNOWN_CONTINUOUS_EVALUATOR, &
    & FHT_UNKNOWN_ENSEMBLE_SOURCE, FHT_UNKNOWN_CONTINUOUS_SOURCE

  PUBLIC FML_INVALID_HANDLE, FML_ERR_NO_ERROR, FML_ERR_UNKNOWN_OBJECT, FML_ERR_INVALID_OBJECT, FML_ERR_INCOMPLETE_OBJECT, &
    & FML_ERR_MISCONFIGURED_OBJECT, FML_ERR_ACCESS_VIOLATION, FML_ERR_FILE_READ, FML_ERR_FILE_WRITE, &
    & FML_ERR_INVALID_PARAMETER_1, FML_ERR_INVALID_PARAMETER_2, FML_ERR_INVALID_PARAMETER_3, FML_ERR_INVALID_PARAMETER_4, &
    & FML_ERR_INVALID_PARAMETER_5, FML_ERR_INVALID_PARAMETER_6, FML_ERR_INVALID_PARAMETER_7, FML_ERR_INVALID_PARAMETER_8, &
    & FML_ERR_UNSUPPORTED

END MODULE FIELDML_API
