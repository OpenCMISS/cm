/* \file
 * $Id: FieldExport.c 542 2009-06-03 17:16:22Z chrispbradley $
 * \author Caton Little
 * \brief 
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is OpenCMISS
 *
 * The Initial Developer of the Original Code is University of Auckland,
 * Auckland, New Zealand and University of Oxford, Oxford, United
 * Kingdom. Portions created by the University of Auckland and University
 * of Oxford are Copyright (C) 2007 by the University of Auckland and
 * the University of Oxford. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */

#ifdef WIN32
#include <windows.h>
#else
#include <stdarg.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FieldExportConstants.h"

/**********************************************************

    Fortran constants. Audit regularly. These aren't in
    FieldExport.h, because the Fortran code #includes
    that file, which leads to a #define collision.

 **********************************************************/
#define BASIS_LINEAR_LAGRANGE_INTERPOLATION     1 //< Linear Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC_LAGRANGE_INTERPOLATION  2 //< Quadratic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_LAGRANGE_INTERPOLATION      3 //< Cubic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_HERMITE_INTERPOLATION       4 //< Cubic Hermite interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC1_HERMITE_INTERPOLATION  5 //< Quadratic Hermite (no derivative at xi=0) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC2_HERMITE_INTERPOLATION  6 //< Quadratic Hermite (no derivative at xi=1) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_LINEAR_SIMPLEX_INTERPOLATION      7 //< Linear Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC_SIMPLEX_INTERPOLATION   8 //< Quadratic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_SIMPLEX_INTERPOLATION       9 //< Cubic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES

#define FIELD_IO_SCALE_FACTORS_NUMBER_TYPE      5
#define FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE    6


#define FIELD_GEOMETRIC_TYPE 1 //Geometric field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_FIBRE_TYPE     2 //Fibre field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_GENERAL_TYPE   3 //General field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_MATERIAL_TYPE  4 //Material field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES

#define FIELD_U_VARIABLE_TYPE    1 //Standard variable type i.e., u \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DELUDELN_VARIABLE_TYPE      2 //Normal derivative variable type i.e., du/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DELUDELT_VARIABLE_TYPE 3 //First time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DEL2UDELT2_VARIABLE_TYPE 4 //Second type derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_V_VARIABLE_TYPE 5 //Second standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DELVDELN_VARIABLE_TYPE 6 //Second normal variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES

#define COORDINATE_RECTANGULAR_CARTESIAN_TYPE 1 //Rectangular Cartesian coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
#define COORDINATE_CYCLINDRICAL_POLAR_TYPE    2 //Cylindrical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
#define COORDINATE_SPHERICAL_POLAR_TYPE       3 //Spherical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
#define COORDINATE_PROLATE_SPHEROIDAL_TYPE    4 //Prolate spheroidal coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
#define COORDINATE_OBLATE_SPHEROIDAL_TYPE     5 //Oblate spheroidal coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES

#define NO_PART_DERIV           1 //No partial derivative i.e., u \see CONSTANTS_PartialDerivativeConstants,CONSTANTS

#define FIRST_PART_DERIV        2 //First partial derivative i.e., du/ds \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define SECOND_PART_DERIV       3 //Second partial derivative i.e., d^2u/ds^2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS

#define PART_DERIV_S1           2 //First partial derivative in the s1 direction i.e., du/ds1 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S1        3 //Second partial derivative in the s1 direction i.e., d^2u/ds1ds1 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S2           4 //First partial derivative in the s2 direction i.e., du/ds2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S2_S2        5 //Second partial derivative in the s2 direction i.e., d^2u/ds2ds2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S2        6 //Cross derivative in the s1 and s2 direction i.e., d^2u/ds1ds2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S3           7 //First partial derivative in the s3 direction i.e., du/ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S3_S3        8 //Second partial derivative in the s3 direction i.e., d^2u/ds3ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S3        9 //Cross derivative in the s1 and s3 direction i.e., d^2u/ds1ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S2_S3       10 //Cross derivative in the s2 and s3 direction i.e., d^2u/ds2ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S2_S3    11 //Cross derivative in the s1, s2 and s3 direction i.e., d^3u/ds1ds2ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S4          12 //First partial derivative in the s4 direction i.e., du/ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S4_S4       13 //Second partial derivative in the s4 direction i.e., d^2u/ds4ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S4       14 //Cross derivative in the s1 and s4 direction i.e., d^2u/ds1ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S2_S4       15 //Cross derivative in the s2 and s4 direction i.e., d^2u/ds2ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S3_S4       16 //Cross derivative in the s3 and s4 direction i.e., d^2u/ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S2_S4    17 //Cross derivative in the s1, s2 and s4 direction i.e., d^3u/ds1ds2ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S3_S4    18 //Cross derivative in the s1, s3 and s4 direction i.e., d^3u/ds1ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S2_S3_S4    19 //Cross derivative in the s2, s3 and s4 direction i.e., d^3u/ds2ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
#define PART_DERIV_S1_S2_S3_S4 20 //Cross derivative in the s1, s2, s3 and s4 direction i.e., d^4u/ds1ds2ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS

/**********************************************************

    Fortran constants. Audit regularly.

 **********************************************************/
typedef struct
{
    int userNumber;
    char isFinished;
    int type;
    int radialInterpolationType;
    int numberOfDimensions;
    double focus;
    double origin[3];
    double orientation[3][3];
}
CMISS_CoordinateSystem;


/*
    API-local structs.
*/
typedef struct
{
    FILE *file;

    int error;
}
FileSession;

typedef struct _SessionListEntry
{
    int handle;

    int type;

    union
    {
        FileSession fileSession;
    };

    struct _SessionListEntry *next;
}
SessionListEntry;

//This is the head of the session list. For convenience, the head itself is a stub, so it exists even for a zero-entry list.
static SessionListEntry sessions;

static int nextHandle = 0;

static int FieldExport_FPrintf( FileSession *const session, const char *format, ... )
{
    va_list args;

    if( session->error == FIELD_EXPORT_NO_ERROR )
    {
        va_start( args, format );
        if( vfprintf( session->file, format, args ) < 0 )
        {
            session->error = FIELD_EXPORT_ERROR_FILE_WRITE;
        }
        va_end( args );
    }

    return session->error;
}


static SessionListEntry *FieldExport_GetSession( const int handle )
{
    SessionListEntry *entry = sessions.next;

    while( entry != NULL )
    {
        if( entry->handle == handle )
        {
            return entry;
        }

        entry = entry->next;
    }

    return NULL;
}


/*
    CMISS-formatted file export routines.
*/
static int FieldExport_File_Group( FileSession *const session, const char *const label )
{
    return FieldExport_FPrintf( session, " Group name: %s\n", label );
}


static int FieldExport_File_MeshDimensions( FileSession *const session, const int dimensions )
{
    return FieldExport_FPrintf( session, " Shape.  Dimension=%d\n", dimensions );
}


static int FieldExport_File_ScalingFactorCount( FileSession *const session, const int scalingFactorCount )
{
    return FieldExport_FPrintf( session, " #Scale factor sets= %d\n", scalingFactorCount );
}


static int FieldExport_File_ScaleFactors( FileSession *const session, const int labelType, const int numberOfXi, const int* const interpolationXi )
{
    int scaleFactorCount = 1;
    int i, j;
    const char * label;
    int *linked = malloc( sizeof(int) * numberOfXi );
    int linkCount = 0;

    for( i = 0; i < numberOfXi; i++ )
    {
        switch( interpolationXi[i] )
        {
        case BASIS_LINEAR_SIMPLEX_INTERPOLATION:
        case BASIS_QUADRATIC_SIMPLEX_INTERPOLATION:
        case BASIS_CUBIC_SIMPLEX_INTERPOLATION:
            linked[i] = 1;
            linkCount++;
            break;
        default:
            linked[i] = 0;
            break;
        }
    }

    FieldExport_FPrintf( session, " " );

    for( i = 0; i < numberOfXi; i++ )
    {
        switch( interpolationXi[i] )
        {
        case BASIS_LINEAR_LAGRANGE_INTERPOLATION:
            scaleFactorCount *= 2;
            label = "l.Lagrange";
            break;
        case BASIS_QUADRATIC_LAGRANGE_INTERPOLATION:
            scaleFactorCount *= 3;
            label = "q.Lagrange";
            break;
        case BASIS_CUBIC_LAGRANGE_INTERPOLATION:
            scaleFactorCount *= 4;
            label = "c.Lagrange";
            break;
        case BASIS_CUBIC_HERMITE_INTERPOLATION:
            scaleFactorCount *= 4;
            label = "c.Hermite";
            break;
        case BASIS_QUADRATIC1_HERMITE_INTERPOLATION:
            scaleFactorCount *= 4;
            label = "q1.Hermite";
            break;
        case BASIS_QUADRATIC2_HERMITE_INTERPOLATION:
            scaleFactorCount *= 4;
            label = "q2.Hermite";
            break;
        case BASIS_LINEAR_SIMPLEX_INTERPOLATION:
            scaleFactorCount = numberOfXi + 1;
            label = "l.simplex";
            break;
        case BASIS_QUADRATIC_SIMPLEX_INTERPOLATION:
            scaleFactorCount = ( numberOfXi + 1 ) * ( numberOfXi + 2 ) / 2;
            label = "q.simplex";
            break;
        case BASIS_CUBIC_SIMPLEX_INTERPOLATION:
            scaleFactorCount = ( numberOfXi + 1 ) * ( numberOfXi + 2 ) * ( numberOfXi + 3 ) / 2;
            label = "c.simplex";
            break;
        default:
            free( linked );
            return FIELD_EXPORT_ERROR_UNKNOWN_INTERPOLATION;
        }

        FieldExport_FPrintf( session, "%s", label );

        if( linkCount > 0 )
        {
            linkCount--;
            FieldExport_FPrintf( session, "(", label );
            for( j = i+1; j < numberOfXi; j++ )
            {
                if( linked[j] == 1 )
                {
                    FieldExport_FPrintf( session, "%d", j+1 );
                }
                if( linkCount > 1 )
                {
                    FieldExport_FPrintf( session, ";" );
                }
                linkCount--;
            }
            FieldExport_FPrintf( session, ")", label );
        }

        if( i < ( numberOfXi - 1 ) )
        {
            FieldExport_FPrintf( session, "*" );
        }
    }
    
    free( linked );

    switch( labelType )
    {
    case FIELD_IO_SCALE_FACTORS_NUMBER_TYPE:
        FieldExport_FPrintf( session, ", #Scale factors=%d\n", scaleFactorCount );
        break;
    case FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE:
        FieldExport_FPrintf( session, ", no modify, standard node based.\n" );
        break;
    default:
        return FIELD_EXPORT_ERROR_UNKNOWN_LABEL_TYPE;
    }

    return session->error;
}


static char *FieldExport_GetVariableLabel( const int fieldType, const int variableType )
{
    switch( fieldType )
    {
    case FIELD_GEOMETRIC_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "unknown";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "Normal_derivative,  field,  normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "first_time_derivative,  field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "second_time_derivative,  field,  second time derivative of variable";
        default:
            return "unknown_geometry,  field,  unknown field variable type";
        }
    case FIELD_FIBRE_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "fibres, anatomical, fibre";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "norm_der_fiber,  normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "first_time_fiber,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "second_time_fiber,  second time derivative of variable";
        default:
            return "unknown_fiber,  unknown field variable type";
        }
    case FIELD_GENERAL_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "general,  field,  rectangular cartesian";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "norm_dev_variable,  field,  rectangular cartesian";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "first_time_variable,  field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "second_time_variable,  field,  second time derivative of variable";
        default:
            return "unknown_general,  field,  unknown field variable type";
        }
    case FIELD_MATERIAL_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "material,  field,  rectangular cartesian";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "normal_material,  field,  normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "fist_time_material,  field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "second_time_material,  field,  second time derivative of variable";
        default:
            return "unknown material,  field,  unknown field variable type";
        }
    default:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "unknown,  field,  unknown standand variable type";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "unknown,  field,  unknown normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "unknown,  field,  unknown first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "unknown, field,  unknown second time derivative of variable";
        default:
            return "unknown,  field,  unknown field variable type";
        }
    }
}


static char *FieldExport_GetCoordinateVariableLabel( CMISS_CoordinateSystem *coordinateSystem )
{
    switch( coordinateSystem->type )
    {
    case COORDINATE_RECTANGULAR_CARTESIAN_TYPE:
        return "coordinates,  coordinate, rectangular cartesian";
    //MUSTDO non-rectangular coordinate systems
    /*
    case COORDINATE_CYCLINDRICAL_POLAR_TYPE:
    case COORDINATE_SPHERICAL_POLAR_TYPE:
    case COORDINATE_PROLATE_SPHEROIDAL_TYPE:
    case COORDINATE_OBLATE_SPHEROIDAL_TYPE:
    */
    default:
        return "unknown";
    }
}


static char *FieldExport_GetCoordinateComponentLabel( CMISS_CoordinateSystem *coordinateSystem, int componentNumber )
{
    switch( coordinateSystem->type )
    {
    case COORDINATE_RECTANGULAR_CARTESIAN_TYPE:
        if( componentNumber == 1 )
        {
            return "x";
        }
        else if( componentNumber == 2 )
        {
            return "y";
        }
        else if( componentNumber == 3 )
        {
            return "z";
        }
        break;
    /*
    case COORDINATE_CYCLINDRICAL_POLAR_TYPE:
    case COORDINATE_SPHERICAL_POLAR_TYPE:
    case COORDINATE_PROLATE_SPHEROIDAL_TYPE:
    case COORDINATE_OBLATE_SPHEROIDAL_TYPE:
    */
    default:
        break;
    }

    return NULL;
}


static int FieldExport_File_NodeCount( FileSession *const session, const int nodeCount )
{
    return FieldExport_FPrintf( session, " #Nodes=           %d\n", nodeCount );
}


static int FieldExport_File_FieldCount( FileSession *const session, const int fieldCount )
{
    return FieldExport_FPrintf( session, " #Fields=%d\n", fieldCount );
}


static int FieldExport_File_CoordinateVariable( FileSession *const session, const int variableIndex,
                                        CMISS_CoordinateSystem *coordinateSystem, const int componentCount )
{
    char *coordinateLabel;

    coordinateLabel = FieldExport_GetCoordinateVariableLabel( coordinateSystem );
    
    return FieldExport_FPrintf( session, " %d) %s, #Components=%d\n", variableIndex, coordinateLabel, componentCount );
}


static int FieldExport_File_Variable( FileSession *const session, const int variableIndex, const int fieldType, const int variableType, const int componentCount )
{
    char *variableLabel;

    variableLabel = FieldExport_GetVariableLabel( fieldType, variableType );
    
    return FieldExport_FPrintf( session, " %d) %s, #Components=%d\n", variableIndex, variableLabel, componentCount );
}


static int FieldExport_File_CoordinateComponent( FileSession *const session, CMISS_CoordinateSystem * coordinateSystem,
    const int componentNumber, const int numberOfXi, const int *const interpolationXi )
{
    const char * const componentLabel = FieldExport_GetCoordinateComponentLabel( coordinateSystem, componentNumber );

    if( componentLabel == NULL )
    {
        FieldExport_FPrintf( session, "   %d.  ", componentNumber );
    }
    else
    {
        FieldExport_FPrintf( session, "   %s.  ", componentLabel );
    }

    if( session->error != FIELD_EXPORT_NO_ERROR )
    {
        return session->error;
    }

    return FieldExport_File_ScaleFactors( session, FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE, numberOfXi, interpolationXi );
}


static int FieldExport_File_Component( FileSession *const session,
    const int componentNumber, const int numberOfXi, const int *const interpolationXi )
{
    if( FieldExport_FPrintf( session, "   %d.  ", componentNumber ) != FIELD_EXPORT_NO_ERROR )
    {
        session->error;
    }

    return FieldExport_File_ScaleFactors( session, FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE, numberOfXi, interpolationXi );
}


static int FieldExport_File_Nodes( FileSession *const session, const int nodeCount, const int *const derivativeCount,
    const int *const elementDerivatives, int firstScaleIndex )
{
    int i, j;
    int derivativeIndex = 0;

    FieldExport_FPrintf( session, "     #Nodes= %d\n", nodeCount );

    for( i = 0; i < nodeCount; i++ )
    {
        FieldExport_FPrintf( session, " %5d.  #Values=%d\n", i + 1, derivativeCount[i] );
        FieldExport_FPrintf( session, "      Value indices:  " );
        for( j = 0; j < derivativeCount[i]; j++ )
        {
            FieldExport_FPrintf( session, " %3d", elementDerivatives[ derivativeIndex++ ] );
        }
        FieldExport_FPrintf( session, "\n" );
        FieldExport_FPrintf( session, "      Scale factor indices: " );
        for( j = 0; j < derivativeCount[i]; j++ )
        {
            firstScaleIndex++;
            FieldExport_FPrintf( session, " %3d", firstScaleIndex );
        }
        FieldExport_FPrintf( session, "\n" );
    }

    return session->error;
}


static int FieldExport_File_ElementIndex( FileSession *session, const int dimensionCount, const int elementIndex )
{
    if( dimensionCount == 3 )
    {
          return FieldExport_FPrintf( session, " Element:            %d 0 0\n", elementIndex );
    }
    else if( dimensionCount == 2 )
    {
          return FieldExport_FPrintf( session, " Element:            0 %d 0\n", elementIndex );
    }
    else
    {
          return FieldExport_FPrintf( session, " Element:            0 0 %d\n", elementIndex );
    }
}


static int FieldExport_File_ElementNodeIndices( FileSession *session, const int nodeCount, const int *const indices )
{
    int i;

    FieldExport_FPrintf( session, "  Nodes:\n" );

    for( i = 0; i < nodeCount; i++ )
    {
        FieldExport_FPrintf( session, "  %10d", indices[i] );
    }

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_File_ElementNodeScales( FileSession *session, const int scaleCount, const double *const scales )
{
    int i;

    FieldExport_FPrintf( session, "   Scale factors:\n" );

    for( i = 0; i < scaleCount; i++ )
    {
        FieldExport_FPrintf( session, "   %.16E", scales[i] );
    }

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_File_OpenSession( const char *const name, int * const handle )
{
    SessionListEntry *session = calloc( 1, sizeof( SessionListEntry ) );

    session->type = EXPORT_TYPE_FILE;
    session->handle = nextHandle++;
    session->fileSession.file = fopen( name, "w" );
    session->fileSession.error = FIELD_EXPORT_NO_ERROR;

    if( session->fileSession.file == NULL )
    {
        free( session );
        return FIELD_EXPORT_ERROR_FILE_IO;
    }

    session->next = sessions.next;
    sessions.next = session;

    *handle = session->handle;

    return FIELD_EXPORT_NO_ERROR;
}


static int FieldExport_File_CloseSession( SessionListEntry *session )
{
    fclose( session->fileSession.file );
    session->type = EXPORT_TYPE_CLOSED;

    return FIELD_EXPORT_NO_ERROR;
}


static int FieldExport_File_NodeNumber( FileSession *session, const int nodeNumber )
{
    return FieldExport_FPrintf( session, " Node:            %d\n", nodeNumber );
}


static int FieldExport_File_NodeValues( FileSession *session, const int valueCount, const double *const values )
{
    int i;

    for( i = 0; i < valueCount; i++ )
    {
        FieldExport_FPrintf( session, "  %.16E", values[i] );
    }
    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_FieldDerivateLabels( FileSession *session, const int numberOfDerivatives, const int *const derivatives )
{
    int i;

    if( ( numberOfDerivatives == 1 ) && ( derivatives[0] == NO_PART_DERIV ) )
    {
        return session->error;
    }

    for( i = 0; i < numberOfDerivatives; i++ )
    {
        switch( derivatives[i] )
        {
        case NO_PART_DERIV:
            break;
        case PART_DERIV_S1:
              FieldExport_FPrintf( session, ", d/ds1" );
              break;
        case PART_DERIV_S1_S1:
              FieldExport_FPrintf( session, ", d2/ds1ds1" );
              break;
        case PART_DERIV_S2:
              FieldExport_FPrintf( session, ", d/ds2" );
              break;
        case PART_DERIV_S2_S2:
              FieldExport_FPrintf( session, ", d2/ds2ds2" );
              break;
        case PART_DERIV_S1_S2:
              FieldExport_FPrintf( session, ", d/ds3" );
              break;
        case PART_DERIV_S3:
              FieldExport_FPrintf( session, ", d2/ds3ds3" );
              break;
        case PART_DERIV_S3_S3:
              FieldExport_FPrintf( session, ", d2/ds3ds3" );
              break;
        case PART_DERIV_S1_S3:
              FieldExport_FPrintf( session, ", d2/ds1ds3" );
              break;
        case PART_DERIV_S2_S3:
              FieldExport_FPrintf( session, ", d2/ds2ds3" );
              break;
        case PART_DERIV_S1_S2_S3:
              FieldExport_FPrintf( session, ", d3/ds1ds2ds3" );
              break;
        case PART_DERIV_S4:
              FieldExport_FPrintf( session, ", d/ds4" );
              break;
        case PART_DERIV_S4_S4:
              FieldExport_FPrintf( session, ", d2/ds4ds4" );
              break;
        case PART_DERIV_S1_S4:
              FieldExport_FPrintf( session, ", d2/ds1ds4" );
              break;
        case PART_DERIV_S2_S4:
              FieldExport_FPrintf( session, ", d2/ds2ds4" );
              break;
        case PART_DERIV_S3_S4:
              FieldExport_FPrintf( session, ", d2/ds3ds4" );
              break;
        case PART_DERIV_S1_S2_S4:
              FieldExport_FPrintf( session, ", d3/ds1ds2ds4" );
              break;
        case PART_DERIV_S1_S3_S4:
              FieldExport_FPrintf( session, ", d3/ds1ds3ds4" );
              break;
        case PART_DERIV_S2_S3_S4:
              FieldExport_FPrintf( session, ", d3/ds2ds3ds4" );
              break;
        case PART_DERIV_S1_S2_S3_S4:
              FieldExport_FPrintf( session, ", d4/ds1ds2ds3ds4" );
              break;
        default:
              FieldExport_FPrintf( session, "unknown field variable type %d", derivatives[i] );
        }
    }

    return session->error;
}


static int FieldExport_File_DerivativeIndices( FileSession *session, const int componentNumber, const int fieldType,
    const int variableType, const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    //MUSTDO add a proper GetComponentLabel( fieldType, variableType, componentNumber ) function.
    FieldExport_FPrintf( session, "   %d.  Value index= %d, #Derivatives= %d", componentNumber, valueIndex, numberOfDerivatives - 1 );

    FieldExport_FieldDerivateLabels( session, numberOfDerivatives, derivatives );

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_File_CoordinateDerivativeIndices( FileSession *session, const int componentNumber, CMISS_CoordinateSystem * coordinateSystem,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    const char * const componentLabel = FieldExport_GetCoordinateComponentLabel( coordinateSystem, componentNumber );

    //MUSTDO add a proper GetComponentLabel( fieldType, variableType, componentNumber ) function.
    if( componentLabel == NULL )
    {
        FieldExport_FPrintf( session, "   %d.  Value index= %d, #Derivatives= %d", componentNumber, valueIndex, numberOfDerivatives - 1 );
    }
    else
    {
        FieldExport_FPrintf( session, "   %s.  Value index= %d, #Derivatives= %d", componentLabel, valueIndex, numberOfDerivatives - 1 );
    }

    FieldExport_FieldDerivateLabels( session, numberOfDerivatives, derivatives );

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


/*
    Public API implementation
*/
int FieldExport_OpenSession( const int type, const char *const name, int * const handle )
{
    if( type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_OpenSession( name, handle );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_Group( const int handle, const char *const label )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_Group( &session->fileSession, label );
    }
    else 
    {
        return FIELD_EXPORT_ERROR_CLOSED_HANDLE;
    }
}


int FieldExport_MeshDimensions( const int handle, const int dimensions )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_MeshDimensions( &session->fileSession, dimensions );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ScalingFactorCount( const int handle, const int scalingFactorCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ScalingFactorCount( &session->fileSession, scalingFactorCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ScaleFactors( const int handle, const int numberOfXi, const int* const interpolationXi )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ScaleFactors( &session->fileSession, FIELD_IO_SCALE_FACTORS_NUMBER_TYPE, numberOfXi, interpolationXi );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_NodeCount( const int handle, const int nodeCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_NodeCount( &session->fileSession, nodeCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_FieldCount( const int handle, const int fieldCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_FieldCount( &session->fileSession, fieldCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_CoordinateVariable( const int handle, const int variableNumber, CMISS_CoordinateSystem coordinateSystem,
    const int componentCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_CoordinateVariable( &session->fileSession, variableNumber, &coordinateSystem, componentCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_Variable( const int handle, const int variableNumber, const int fieldType, const int variableType,
    const int componentCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_Variable( &session->fileSession, variableNumber, fieldType, variableType, componentCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_CoordinateComponent( const int handle, CMISS_CoordinateSystem coordinateSystem,
    const int componentNumber, const int numberOfXi, const int * const interpolationXi )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_CoordinateComponent( &session->fileSession, &coordinateSystem, componentNumber, numberOfXi, interpolationXi );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_Component( const int handle, const int componentNumber, const int numberOfXi, const int * const interpolationXi )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_Component( &session->fileSession, componentNumber, numberOfXi, interpolationXi );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_Nodes( const int handle, const int nodeCount, const int *const derivativeCount,
    const int *const elementDerivatives, const int firstScaleIndex )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_Nodes( &session->fileSession, nodeCount, derivativeCount, elementDerivatives, firstScaleIndex );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ElementIndex( const int handle, const int dimensionCount, const int index )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ElementIndex( &session->fileSession, dimensionCount, index );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ElementNodeIndices( const int handle, const int nodeCount, const int* const indices )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ElementNodeIndices( &session->fileSession, nodeCount, indices );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ElementNodeScales( const int handle, const int scaleCount, const double* const scales )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ElementNodeScales( &session->fileSession, scaleCount, scales );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_CloseSession( const int handle )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_CloseSession( session );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}


int FieldExport_NodeNumber( const int handle, const int nodeNumber )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_NodeNumber( &session->fileSession, nodeNumber );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}


int FieldExport_NodeValues( const int handle, const int valueCount, const double* const values )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_NodeValues( &session->fileSession, valueCount, values );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}


int FieldExport_CoordinateDerivativeIndices( const int handle, const int componentNumber, CMISS_CoordinateSystem coordinateSystem,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_CoordinateDerivativeIndices( &session->fileSession, componentNumber, &coordinateSystem, numberOfDerivatives, derivatives, valueIndex );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}


int FieldExport_DerivativeIndices( const int handle, const int componentNumber, const int fieldType, const int variableType,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_DerivativeIndices( &session->fileSession, componentNumber, fieldType, variableType, numberOfDerivatives, derivatives, valueIndex );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}
