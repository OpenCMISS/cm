/* \file
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
 * Auckland, New Zealand, the University of Oxford, Oxford, United
 * Kingdom and King's College, London, United Kingdom. Portions created
 * by the University of Auckland, the University of Oxford and King's
 * College, London are Copyright (C) 2007-2010 by the University of
 * Auckland, the University of Oxford and King's College, London.
 * All Rights Reserved.
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

//#include <hdf5.h>
#ifdef H5_VERS_MAJOR
#define USE_HDF5
#endif

#include "FieldExportConstants.h"

/**********************************************************

    Fortran constants. Audit regularly. These aren't in
    FieldExport.h, because the Fortran code #includes
    that file, which leads to a #define collision.

 **********************************************************/

#define BASIS_LAGRANGE_HERMITE_TP_TYPE          1 //< Lagrange Hermite tensor product basis \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
#define BASIS_SIMPLEX_TYPE                      2 //< Simplex basis \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES


#define BASIS_LINEAR_LAGRANGE_INTERPOLATION     1 //< Linear Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC_LAGRANGE_INTERPOLATION  2 //< Quadratic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_LAGRANGE_INTERPOLATION      3 //< Cubic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_HERMITE_INTERPOLATION       4 //< Cubic Hermite interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC1_HERMITE_INTERPOLATION  5 //< Quadratic Hermite (no derivative at xi=0) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC2_HERMITE_INTERPOLATION  6 //< Quadratic Hermite (no derivative at xi=1) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_LINEAR_SIMPLEX_INTERPOLATION      7 //< Linear Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC_SIMPLEX_INTERPOLATION   8 //< Quadratic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_SIMPLEX_INTERPOLATION       9 //< Cubic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES

#define FIELD_IO_INTERPOLATION_HEADER_SCALE     4
#define FIELD_IO_INTERPOLATION_HEADER_NODAL     5
#define FIELD_IO_INTERPOLATION_HEADER_GRID      6
#define FIELD_IO_INTERPOLATION_HEADER_GAUSS     7
#define FIELD_IO_INTERPOLATION_HEADER_CONSTANT  8


#define FIELD_GEOMETRIC_TYPE 1 //Geometric field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_FIBRE_TYPE     2 //Fibre field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_GENERAL_TYPE   3 //General field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_MATERIAL_TYPE  4 //Material field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_GEOMETRIC_GENERAL_TYPE 5 //Geometric general field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES

#define FIELD_U_VARIABLE_TYPE    1 //Standard variable type i.e., u \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DELUDELN_VARIABLE_TYPE      2 //Normal derivative variable type i.e., du/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DELUDELT_VARIABLE_TYPE 3 //First time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DEL2UDELT2_VARIABLE_TYPE 4 //Second type derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_V_VARIABLE_TYPE 5 //Second standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_DELVDELN_VARIABLE_TYPE 6 //Second normal variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_U1_VARIABLE_TYPE 9 //Third standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_U2_VARIABLE_TYPE 13 //Fourth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES

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

/*
    API-local structs.
*/
typedef struct
{
    FILE *file;
#ifdef USE_HDF5
    hid_t hd5Handle;
#endif

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

#ifdef USE_HDF5
static void eep()
{
    int *i = NULL;
    int j;

    FILE *f = fopen( "D:\\err.txt", "w" );
    H5Eprint1(f);
    fclose(f);

    j = *i;
}
#endif


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

int FieldExport_InterpolationType( const int interpType );

/*
    CMISS-formatted file export routines.
*/
static int FieldExport_File_Group( FileSession *const session, const char *const label )
{
  return FieldExport_FPrintf( session, " Group name: %s\n", label );
    //return FieldExport_FPrintf( session, " Region: /%s\n", label );
}


static int FieldExport_File_MeshDimensions( FileSession *const session, const int dimensions, const int basisType )
{
  switch( basisType )
    {
    case BASIS_LAGRANGE_HERMITE_TP_TYPE:
      return FieldExport_FPrintf( session, " Shape.  Dimension=%d\n", dimensions );
    case BASIS_SIMPLEX_TYPE:
      switch( dimensions )
	{
	case 1:
	  return FieldExport_FPrintf( session, " Shape.  Dimension=%d\n", dimensions );
	case 2:
	  return FieldExport_FPrintf( session, " Shape.  Dimension=%d, simplex(2)*simplex\n", dimensions );
	case 3:
	  return FieldExport_FPrintf( session, " Shape.  Dimension=%d, simplex(2;3)*simplex*simplex\n", dimensions );
	default:
	  return FieldExport_FPrintf( session, " Shape.  Dimension=%d\n", dimensions );
	}
    default:
      return FieldExport_FPrintf( session, " Shape.  Dimension=%d\n", dimensions );
    }
}

static int FieldExport_File_ScalingFactorCount( FileSession *const session, const int scalingFactorCount )
{
    return FieldExport_FPrintf( session, " #Scale factor sets= %d\n", scalingFactorCount );
}


static int FieldExport_File_InterpolationHeader( FileSession *const session, const int labelType, const int numberOfXi, const int* const interpolationXi)
{
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
        if( labelType == FIELD_IO_INTERPOLATION_HEADER_GRID )
        {
	  label = "l.Lagrange";
        }
        else if( labelType == FIELD_IO_INTERPOLATION_HEADER_CONSTANT )
        {
            label = "constant";
        }
        else if( labelType == FIELD_IO_INTERPOLATION_HEADER_GAUSS )
        {
	  label = "l.Lagrange"; 
        }
        else
        {
            switch( interpolationXi[i] )
            {
            case BASIS_LINEAR_LAGRANGE_INTERPOLATION:
                label = "l.Lagrange";
                break;
            case BASIS_QUADRATIC_LAGRANGE_INTERPOLATION:
                label = "q.Lagrange";
                break;
            case BASIS_CUBIC_LAGRANGE_INTERPOLATION:
                label = "c.Lagrange";
                break;
            case BASIS_CUBIC_HERMITE_INTERPOLATION:
                label = "c.Hermite";
                break;
            case BASIS_QUADRATIC1_HERMITE_INTERPOLATION:
                label = "LagrangeHermite";
                break;
            case BASIS_QUADRATIC2_HERMITE_INTERPOLATION:
                label = "HermiteLagrange";
                break;
            case BASIS_LINEAR_SIMPLEX_INTERPOLATION:
                label = "l.simplex";
                break;
            case BASIS_QUADRATIC_SIMPLEX_INTERPOLATION:
                label = "q.simplex";
                break;
            case BASIS_CUBIC_SIMPLEX_INTERPOLATION:
                label = "c.simplex";
                break;
            default:
                free( linked );
                return FIELD_EXPORT_ERROR_UNKNOWN_INTERPOLATION;
            }
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
    case FIELD_IO_INTERPOLATION_HEADER_NODAL:
        FieldExport_FPrintf( session, ", no modify, standard node based.\n" );
        break;
    case FIELD_IO_INTERPOLATION_HEADER_GRID:
        FieldExport_FPrintf( session, ", no modify, grid based.\n" );
        break;
    case FIELD_IO_INTERPOLATION_HEADER_CONSTANT:
        FieldExport_FPrintf( session, ", no modify, grid based.\n" );
        break;
    case FIELD_IO_INTERPOLATION_HEADER_GAUSS:
        FieldExport_FPrintf( session, ", no modify, grid based.\n" );
        break;
    default:
        return FIELD_EXPORT_ERROR_UNKNOWN_LABEL_TYPE;
    }

    return session->error;
}


static int FieldExport_File_InterpolationHeaderScale( FileSession *const session, const int numberOfXi, const int* const interpolationXi,
        const int numberOfScaleFactors)
{
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
            label = "l.Lagrange";
            break;
        case BASIS_QUADRATIC_LAGRANGE_INTERPOLATION:
            label = "q.Lagrange";
            break;
        case BASIS_CUBIC_LAGRANGE_INTERPOLATION:
            label = "c.Lagrange";
            break;
        case BASIS_CUBIC_HERMITE_INTERPOLATION:
            label = "c.Hermite";
            break;
        case BASIS_QUADRATIC1_HERMITE_INTERPOLATION:
            label = "LagrangeHermite";
            break;
        case BASIS_QUADRATIC2_HERMITE_INTERPOLATION:
            label = "HermiteLagrange";
            break;
        case BASIS_LINEAR_SIMPLEX_INTERPOLATION:
            label = "l.simplex";
            break;
        case BASIS_QUADRATIC_SIMPLEX_INTERPOLATION:
            label = "q.simplex";
            break;
        case BASIS_CUBIC_SIMPLEX_INTERPOLATION:
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

    FieldExport_FPrintf( session, ", #Scale factors=%d\n", numberOfScaleFactors );

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
            return "field,  normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "field,  second time derivative of variable";
        default:
            return "field,  real";
        }
    case FIELD_FIBRE_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "anatomical, fibre";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "second time derivative of variable";
        default:
            return "real";
        }
    case FIELD_GENERAL_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_V_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_U1_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_U2_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "field,  second time derivative of variable";
        default:
            return "field,  real";
        }
    case FIELD_MATERIAL_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "field,  normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "field,  second time derivative of variable";
        default:
            return "field,  real";
        }
    case FIELD_GEOMETRIC_GENERAL_TYPE:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "field,  rectangular cartesian";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "field,  first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "field,  second time derivative of variable";
        default:
            return "field,  real";
        }
    default:
        switch( variableType )
        {
        case FIELD_U_VARIABLE_TYPE:
            return "field,  unknown standand variable type";
        case FIELD_DELUDELN_VARIABLE_TYPE:
            return "field,  unknown normal derivative of variable";
        case FIELD_DELUDELT_VARIABLE_TYPE:
            return "field,  unknown first time derivative of variable";
        case FIELD_DEL2UDELT2_VARIABLE_TYPE:
            return "field,  unknown second time derivative of variable";
        default:
            return "field,  real";
        }
    }
}


static char *FieldExport_GetCoordinateVariableLabel( int coordinateSystemType )
{
    switch( coordinateSystemType )
    {
    case COORDINATE_RECTANGULAR_CARTESIAN_TYPE:
        return "coordinate, rectangular cartesian";
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


static char *FieldExport_GetCoordinateComponentLabel( int coordinateSystemType, int componentNumber )
{
    switch( coordinateSystemType )
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


static int FieldExport_File_CoordinateVariable( FileSession *const session, const char *variableName, const int variableIndex,
                                        int coordinateSystemType, const int componentCount )
{
    char *coordinateLabel;

    coordinateLabel = FieldExport_GetCoordinateVariableLabel( coordinateSystemType );
    
    return FieldExport_FPrintf( session, " %d) %s, %s, #Components=%d\n", variableIndex, variableName, coordinateLabel, componentCount );
}


static int FieldExport_File_Variable( FileSession *const session, const char *variableName, const int variableIndex, 
				      const int fieldType, const int variableType, const int componentCount )
{
    char *variableLabel;

    variableLabel = FieldExport_GetVariableLabel( fieldType, variableType );
    
    return FieldExport_FPrintf( session, " %d) %s, %s, #Components=%d\n", variableIndex, variableName, variableLabel, componentCount );
}


static int FieldExport_File_CoordinateComponent( FileSession *const session, int coordinateSystemType,
    const int componentNumber, const int interpType, const int numberOfXi, const int *const interpolationXi )
{
    const char * const componentLabel = FieldExport_GetCoordinateComponentLabel( coordinateSystemType, componentNumber );
    int headerType;

    headerType = FieldExport_InterpolationType( interpType );
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

    return FieldExport_File_InterpolationHeader( session, headerType, numberOfXi, interpolationXi );
}


static int FieldExport_File_Component( FileSession *const session,
    const int componentNumber, const int interpType, const int numberOfXi, const int *const interpolationXi )
{
  const int headerType = FieldExport_InterpolationType( interpType );

    if( FieldExport_FPrintf( session, "   %d.  ", componentNumber ) != FIELD_EXPORT_NO_ERROR )
    {
        return session->error;
    }

    return FieldExport_File_InterpolationHeader( session, headerType, numberOfXi, interpolationXi );
}


static int FieldExport_File_ElementGridSize( FileSession *const session, const int interpType, const int numberOfXi, const int *const numberGauss )
{
  int i,numGrid[3];
  const int headerType = FieldExport_InterpolationType( interpType );

  if( headerType == FIELD_IO_INTERPOLATION_HEADER_CONSTANT)
    {
      numGrid[0]=0;
      numGrid[1]=0;
      numGrid[2]=0;
    }
  else if( headerType == FIELD_IO_INTERPOLATION_HEADER_GRID )
    {
      numGrid[0]=1;
      numGrid[1]=1;
      numGrid[2]=1;
    }
  else if( headerType == FIELD_IO_INTERPOLATION_HEADER_GAUSS )
    {
    for( i = 0; i < numberOfXi; i++ )
      {
	numGrid[i]=numberGauss[i]-1;
      }
    }
  else
    {
      numGrid[0]=1;
      numGrid[1]=1;
      numGrid[2]=1;
    }

    if( FieldExport_FPrintf( session, "     " ) != FIELD_EXPORT_NO_ERROR )
    {
        return session->error;
    }

    for( i = 0; i < numberOfXi; i++ )
    {
      if( FieldExport_FPrintf( session, "#xi%d=%d", i+1, numGrid[i] ) != FIELD_EXPORT_NO_ERROR )
        {
            return session->error;
        }

        if( i < ( numberOfXi - 1 ) )
        {
            if( FieldExport_FPrintf( session, ", " ) != FIELD_EXPORT_NO_ERROR )
            {
                return session->error;
            }
        }
    }
    if( FieldExport_FPrintf( session, "\n" ) != FIELD_EXPORT_NO_ERROR )
    {
        return session->error;
    }

    return session->error;
}


static int FieldExport_File_NodeScaleIndexes( FileSession *const session, const int nodeCount, const int *const derivativeCount,
    const int *const elementDerivatives, const int *const nodeIndexes, const int *const scaleIndexes )
{
    int i, j;
    int derivativeIndex = 0;
    int scaleIndex = 0;

    FieldExport_FPrintf( session, "     #Nodes= %d\n", nodeCount );

    for( i = 0; i < nodeCount; i++ )
    {
        FieldExport_FPrintf( session, " %5d.  #Values=%d\n", nodeIndexes[i], derivativeCount[i] );
        FieldExport_FPrintf( session, "      Value indices:  " );
        for( j = 0; j < derivativeCount[i]; j++ )
        {
            FieldExport_FPrintf( session, " %3d", elementDerivatives[ derivativeIndex++ ] );
        }
        FieldExport_FPrintf( session, "\n" );
        FieldExport_FPrintf( session, "      Scale factor indices: " );
        for( j = 0; j < derivativeCount[i]; j++ )
        {
            //We're currently using firstScaleIndex == -1 to tell us that there's no scaling. Ugly but functional.
            if( scaleIndex >= 0 )
            {
                FieldExport_FPrintf( session, " %3d", scaleIndexes[scaleIndex++] );
            }
            else
            {
                FieldExport_FPrintf( session, " 0" );
            }
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
      /* return FieldExport_FPrintf( session, " Element:            0 %d 0\n", elementIndex ); */
          return FieldExport_FPrintf( session, " Element:            %d 0 0\n", elementIndex );
    }
    else
    {
      /* return FieldExport_FPrintf( session, " Element:            0 0 %d\n", elementIndex ); */
       return FieldExport_FPrintf( session, " Element:            %d 0 0\n", elementIndex ); 
    }
}


static int FieldExport_File_ElementNodeIndices( FileSession *session, const int nodeCount, const int *const indices )
{
    int i;

    FieldExport_FPrintf( session, "   Nodes:\n" );

    for( i = 0; i < nodeCount; i++ )
    {
        FieldExport_FPrintf( session, "  %10d", indices[i] );
    }

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_File_ElementNodeScales( FileSession *session, const int isFirstSet, const int scaleCount, const double *const scales )
{
    int i;

    if( isFirstSet )
    {
        FieldExport_FPrintf( session, "   Scale factors:\n" );
    }

    FieldExport_FPrintf( session, "  " );

    for( i = 0; i < scaleCount; i++ )
    {
        FieldExport_FPrintf( session, "   %.16E", scales[i] );
    }

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_File_ElementGridValues( FileSession *session, const int isFirstSet, const int valueCount, const double value )
{
    int i;
    /* int valueCount = 1 << dimensionCount; */

    if( isFirstSet )
    {
        FieldExport_FPrintf( session, "   Values:\n" );
    }

    FieldExport_FPrintf( session, "  " );

    for( i = 0; i < valueCount; i++ )
    {
        FieldExport_FPrintf( session, "   %.16E", value );
    }

    FieldExport_FPrintf( session, "\n" );

    return session->error;
}


static int FieldExport_File_OpenSession( const char *const name, int * const handle )
{
    SessionListEntry *session = calloc( 1, sizeof( SessionListEntry ) );
    char hd5Name[256];

    strcpy( hd5Name, name );
    strcat( hd5Name, ".h5" );

    session->type = EXPORT_TYPE_FILE;
    session->handle = nextHandle++;
    session->fileSession.file = fopen( name, "w" );

#ifdef USE_HDF5
    session->fileSession.hd5Handle = H5Fcreate( hd5Name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
#endif

    session->fileSession.error = FIELD_EXPORT_NO_ERROR;

    if( ( session->fileSession.file == NULL )
#ifdef USE_HDF5
        || ( session->fileSession.hd5Handle < 0 )
#endif
        )
    {
        if( session->fileSession.file != NULL )
        {
            fclose( session->fileSession.file );
        }
#ifdef USE_HDF5
        if( session->fileSession.hd5Handle >= 0 )
        {
            H5Fclose( session->fileSession.hd5Handle );
        }
#endif
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

#ifdef USE_HDF5
    H5Fclose( session->fileSession.hd5Handle );
#endif

    session->type = EXPORT_TYPE_CLOSED;

    return FIELD_EXPORT_NO_ERROR;
}


#ifdef USE_HDF5
static int FieldExport_File_HD5_NodeValues( FileSession *session, const int nodeNumber, const int valueCount, const double *const values )
{
    hid_t dataset_id, dataspace_id, attribute_id, attribute_dataspace_id;
    hsize_t dims[1];
    int attributes[1];
    herr_t status;
    char setName[256];

    dims[0] = valueCount; 
    if( ( dataspace_id = H5Screate_simple( 1, dims, NULL ) ) < 0 )
    {
        return dataspace_id;
    }

    sprintf( setName, "/node%d", nodeNumber );
    if( ( dataset_id = H5Dcreate( session->hd5Handle, setName, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) ) < 0 )
    {
        return dataset_id;
    }

    dims[0] = 1; 
    if( ( attribute_dataspace_id = H5Screate_simple( 1, dims, NULL ) ) < 0 )
    {
        return dataspace_id;
    }
    if( ( attribute_id = H5Acreate( dataset_id, "Node", H5T_STD_I32BE, attribute_dataspace_id, H5P_DEFAULT, H5P_DEFAULT ) ) < 0 )
    {
        return attribute_id;
    }

    attributes[0] = nodeNumber;
    if( ( status = H5Awrite( attribute_id, H5T_NATIVE_INT, attributes ) ) < 0 )
    {
        return status;
    }

    if( ( status = H5Aclose( attribute_id ) ) < 0 )
    {
        return status;
    }

    if( ( status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values ) ) < 0 )
    {
        return status;
    }

    if( ( status = H5Dclose( dataset_id ) ) < 0 )
    {
        return status;
    }

    if( ( status = H5Sclose( attribute_dataspace_id ) ) < 0 )
    {
        return status;
    }

    if( ( status = H5Sclose( dataspace_id ) ) < 0 )
    {
        return status;
    }

    return 0;
}
#endif

static int FieldExport_File_NodeValues( FileSession *session, const int nodeNumber, const int valueCount, const double *const values )
{
    int i;
#ifdef USE_HDF5
    herr_t status;
#endif
    static int lastNodeNumber = -1; //A little bit of a hack, but then so is the whole file format.

    if( nodeNumber != lastNodeNumber )
    {
        lastNodeNumber = nodeNumber;
        FieldExport_FPrintf( session, " Node:            %d\n", nodeNumber );
    }

    for( i = 0; i < valueCount; i++ )
    {
        FieldExport_FPrintf( session, "  %.16E", values[i] );
    }
    FieldExport_FPrintf( session, "\n" );

#ifdef USE_HDF5
    status = FieldExport_File_HD5_NodeValues( session, nodeNumber, valueCount, values );

    if( status < 0 )
    {
        eep();
        session->error = FIELD_EXPORT_ERROR_HDF5_ERROR;
    }
#endif

    return session->error;
}


static int FieldExport_FieldDerivateLabels( FileSession *session, const int numberOfDerivatives, const int *const derivatives )
{
    int i;

    if( ( numberOfDerivatives == 1 ) && ( derivatives[0] == NO_PART_DERIV ) )
    {
        return session->error;
    }

    FieldExport_FPrintf( session, "(" );
    for( i = 0; i < numberOfDerivatives; i++ )
    {
        if( i > 1 )
        {
            FieldExport_FPrintf( session, "," );
        }
        switch( derivatives[i] )
        {
        case NO_PART_DERIV:
            break;
        case PART_DERIV_S1:
              FieldExport_FPrintf( session, "d/ds1" );
              break;
        case PART_DERIV_S1_S1:
              FieldExport_FPrintf( session, "d2/ds1ds1" );
              break;
        case PART_DERIV_S2:
              FieldExport_FPrintf( session, "d/ds2" );
              break;
        case PART_DERIV_S2_S2:
              FieldExport_FPrintf( session, "d2/ds2ds2" );
              break;
        case PART_DERIV_S3:
              FieldExport_FPrintf( session, "d/ds3" );
              break;
        case PART_DERIV_S3_S3:
              FieldExport_FPrintf( session, "d2/ds3ds3" );
              break;
        case PART_DERIV_S1_S2:
              FieldExport_FPrintf( session, "d2/ds1ds2" );
              break;
        case PART_DERIV_S1_S3:
              FieldExport_FPrintf( session, "d2/ds1ds3" );
              break;
        case PART_DERIV_S2_S3:
              FieldExport_FPrintf( session, "d2/ds2ds3" );
              break;
        case PART_DERIV_S1_S2_S3:
              FieldExport_FPrintf( session, "d3/ds1ds2ds3" );
              break;
        case PART_DERIV_S4:
              FieldExport_FPrintf( session, "d/ds4" );
              break;
        case PART_DERIV_S4_S4:
              FieldExport_FPrintf( session, "d2/ds4ds4" );
              break;
        case PART_DERIV_S1_S4:
              FieldExport_FPrintf( session, "d2/ds1ds4" );
              break;
        case PART_DERIV_S2_S4:
              FieldExport_FPrintf( session, "d2/ds2ds4" );
              break;
        case PART_DERIV_S3_S4:
              FieldExport_FPrintf( session, "d2/ds3ds4" );
              break;
        case PART_DERIV_S1_S2_S4:
              FieldExport_FPrintf( session, "d3/ds1ds2ds4" );
              break;
        case PART_DERIV_S1_S3_S4:
              FieldExport_FPrintf( session, "d3/ds1ds3ds4" );
              break;
        case PART_DERIV_S2_S3_S4:
              FieldExport_FPrintf( session, "d3/ds2ds3ds4" );
              break;
        case PART_DERIV_S1_S2_S3_S4:
              FieldExport_FPrintf( session, "d4/ds1ds2ds3ds4" );
              break;
        default:
              FieldExport_FPrintf( session, "real" );
        }
    }
    FieldExport_FPrintf( session, ")" );

    return session->error;
}


static int FieldExport_File_DerivativeIndices( FileSession *session, const int componentNumber, const int fieldType,
    const int variableType, const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    //MUSTDO add a proper GetComponentLabel( fieldType, variableType, componentNumber ) function.
    FieldExport_FPrintf( session, "   %d.  Value index= %d, #Derivatives= %d", componentNumber, valueIndex, numberOfDerivatives - 1 );

    FieldExport_FieldDerivateLabels( session, numberOfDerivatives, derivatives );

    return session->error;
}


static int FieldExport_File_CoordinateDerivativeIndices( FileSession *session, const int componentNumber, int coordinateSystemType,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    const char * const componentLabel = FieldExport_GetCoordinateComponentLabel( coordinateSystemType, componentNumber );

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

    return session->error;
}


/*
    Public API implementation
*/

int FieldExport_InterpolationType( const int interpType )
{
    if(interpType == 1) 
      {
	return  FIELD_IO_INTERPOLATION_HEADER_CONSTANT;
      }
    else if(interpType == 2)
      {
	return FIELD_IO_INTERPOLATION_HEADER_CONSTANT;
      }
    else if(interpType == 3)
      {
	return FIELD_IO_INTERPOLATION_HEADER_NODAL;
      }
    else if(interpType == 4)
      {
	return FIELD_IO_INTERPOLATION_HEADER_GRID;
      }
    else if(interpType == 5)
      {
	return FIELD_IO_INTERPOLATION_HEADER_GAUSS;
      }
    else if(interpType == 6)
      {
	return FIELD_IO_INTERPOLATION_HEADER_NODAL;
      }
    else
      {
	return FIELD_IO_INTERPOLATION_HEADER_GRID;
      }

}


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


int FieldExport_MeshDimensions( const int handle, const int dimensions, const int basisType )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
      return FieldExport_File_MeshDimensions( &session->fileSession, dimensions, basisType );
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


int FieldExport_ScaleFactors( const int handle, const int numberOfXi, const int* const interpolationXi, const int numberOfScaleFactors )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_InterpolationHeaderScale( &session->fileSession, numberOfXi, interpolationXi, numberOfScaleFactors );
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


int FieldExport_CoordinateVariable( const int handle, const char *variableName, const int variableNumber, int coordinateSystemType,
    const int componentCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
      return FieldExport_File_CoordinateVariable( &session->fileSession, variableName, variableNumber, coordinateSystemType, componentCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_Variable( const int handle, const char *variableName, const int variableNumber, const int fieldType, const int variableType,
    const int componentCount )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
      return FieldExport_File_Variable( &session->fileSession, variableName, variableNumber, fieldType, variableType, componentCount );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_CoordinateComponent( const int handle, int coordinateSystemType,
    const int componentNumber, const int interpType, const int numberOfXi, const int * const interpolationXi )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_CoordinateComponent( &session->fileSession, coordinateSystemType, componentNumber, interpType, numberOfXi, interpolationXi );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_Component( const int handle, const int componentNumber, const int interpType, const int numberOfXi, const int * const interpolationXi )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_Component( &session->fileSession, componentNumber, interpType, numberOfXi, interpolationXi );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ElementGridSize( const int handle, const int interpType, const int numberOfXi, const int *const numberGauss  )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
      return FieldExport_File_ElementGridSize( &session->fileSession, interpType, numberOfXi, numberGauss );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_NodeScaleIndexes( const int handle, const int nodeCount, const int *const derivativeCount,
    const int *const elementDerivatives, const int *const nodeIndexes, const int *const scaleIndexes )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_NodeScaleIndexes( &session->fileSession, nodeCount, derivativeCount, elementDerivatives, nodeIndexes, scaleIndexes );
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


int FieldExport_ElementNodeScales( const int handle, const int isFirstSet, const int scaleCount, const double* const scales )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ElementNodeScales( &session->fileSession, isFirstSet, scaleCount, scales );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }
}


int FieldExport_ElementGridValues( const int handle, const int isFirstSet, const int dimensionCount, const double value )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_ElementGridValues( &session->fileSession, isFirstSet, dimensionCount, value );
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


int FieldExport_NodeValues( const int handle, const int nodeNumber, const int valueCount, const double* const values )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_NodeValues( &session->fileSession, nodeNumber, valueCount, values );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}


int FieldExport_CoordinateDerivativeIndices( const int handle, const int componentNumber, const int coordinateSystemType,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
    SessionListEntry *session = FieldExport_GetSession( handle );
    
    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_File_CoordinateDerivativeIndices( &session->fileSession, componentNumber, coordinateSystemType, numberOfDerivatives, derivatives, valueIndex );
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

int FieldExport_VersionInfo( const int handle, const int numberOfVersions )
{
    SessionListEntry *session = FieldExport_GetSession( handle );

    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        if (numberOfVersions > 1) {
            return FieldExport_FPrintf( &session->fileSession, ", #Versions=%d", numberOfVersions );
        }
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    return FIELD_EXPORT_NO_ERROR;
}

int FieldExport_EndComponent( const int handle )
{
    SessionListEntry *session = FieldExport_GetSession( handle );

    if( session == NULL )
    {
        return FIELD_EXPORT_ERROR_BAD_HANDLE;
    }
    else if( session->type == EXPORT_TYPE_FILE )
    {
        return FieldExport_FPrintf( &session->fileSession, "\n" );
    }
    else
    {
        return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
    }

    /* Shouldn't get to here */
    return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
}
