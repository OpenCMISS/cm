#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FieldExportConstants.h"

/**********************************************************

	Fortran constants. Audit regularly. These aren't in
	FieldExport.h, because the Fortran code #includes
	that file, which leads to a #define collision.

 **********************************************************/
#define BASIS_LINEAR_LAGRANGE_INTERPOLATION		1 //< Linear Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC_LAGRANGE_INTERPOLATION	2 //< Quadratic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_LAGRANGE_INTERPOLATION		3 //< Cubic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_HERMITE_INTERPOLATION		4 //< Cubic Hermite interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC1_HERMITE_INTERPOLATION	5 //< Quadratic Hermite (no derivative at xi=0) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC2_HERMITE_INTERPOLATION	6 //< Quadratic Hermite (no derivative at xi=1) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_LINEAR_SIMPLEX_INTERPOLATION		7 //< Linear Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_QUADRATIC_SIMPLEX_INTERPOLATION	8 //< Quadratic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
#define BASIS_CUBIC_SIMPLEX_INTERPOLATION		9 //< Cubic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES

#define FIELD_IO_SCALE_FACTORS_NUMBER_TYPE		5
#define FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE	6


#define FIELD_GEOMETRIC_TYPE 1 //Geometric field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_FIBRE_TYPE     2 //Fibre field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_GENERAL_TYPE   3 //General field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
#define FIELD_MATERIAL_TYPE  4 //Material field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES

#define FIELD_STANDARD_VARIABLE_TYPE    1 //Standard variable type i.e., u \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_NORMAL_VARIABLE_TYPE      2 //Normal derivative variable type i.e., du/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_TIME_DERIV1_VARIABLE_TYPE 3 //First time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
#define FIELD_TIME_DERIV2_VARIABLE_TYPE 4 //Second type derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES

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


#define HANDLE_TYPE_CHUNK_SIZE 32

/*
	API-local structs.
*/
typedef struct
{
	int handle;

	FILE *file;
}
FileSession;

typedef struct _FileSessionListEntry
{
	FileSession *session;

	struct _FileSessionListEntry *next;
}
FileSessionListEntry;

static FileSessionListEntry fileSessions;

static int nextHandle = 0;

static int *handleTypes = NULL;


static void flog( char *logline, int line )
{
	char foo[1024];
	FILE *f;

	sprintf( foo, "c:\\log%d.bob", GetCurrentThreadId() );
	f = fopen( foo, "a" );
	fprintf( f, "%d: %s\n", line, logline );
	fclose(f);
}


static void fblog( char *string, int data )
{
	char foo[1024];
	FILE *f;

	sprintf( foo, "c:\\log%d.bob", GetCurrentThreadId() );
	f = fopen( foo, "a" );
	fprintf( f, "%s %08x\n", string, data );
	fclose(f);
}


/*
	General-purpose handle routines
*/
static int FieldExport_GetNewSessionHandle( const int type )
{
	if( ( nextHandle % HANDLE_TYPE_CHUNK_SIZE ) == 0 )
	{
		int *newHandleTypes = calloc( nextHandle + HANDLE_TYPE_CHUNK_SIZE, sizeof( int ) );

		if( handleTypes != NULL )
		{
			memcpy( newHandleTypes, handleTypes, nextHandle * sizeof( int ) );
		}

		handleTypes = newHandleTypes;
	}

	handleTypes[nextHandle] = type;

	return nextHandle++;
}


static int FieldExport_CheckHandle( const int handle )
{
	if( ( handle < 0 ) || ( handle >= nextHandle ) )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}
	else if( handleTypes[handle] == EXPORT_TYPE_CLOSED )
	{
		return FIELD_EXPORT_ERROR_CLOSED_HANDLE;
	}

	return FIELD_EXPORT_NO_ERROR;
}


/*
	CMISS-formatted file export routines.
*/
static int FieldExport_File_Group( const FileSession *const session, const char *const label )
{
	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	if( fprintf( session->file, " Group name: %s\n", label ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_MeshDimensions( const FileSession *const session, const int dimensions )
{
	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	if( fprintf( session->file, " Shape.  Dimension=%d\n", dimensions ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_ScalingFactorCount( const FileSession *const session, const int scalingFactorCount )
{
	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	if( fprintf( session->file, " #Scale factor sets= %d\n", scalingFactorCount ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_LagrangeHermiteScaleFactors( const FileSession *const session, const int labelType, const int numberOfXi, const int* const interpolationXi )
{
	int err;
	int scaleFactorCount = 1;
	int i;

	flog( __FILE__, __LINE__ );

	fblog( "interpolationXi = ", interpolationXi );

	fflush(NULL);

	fprintf( session->file, " " );

	//MUSTDO Extract this loop, and write to a string.
	for( i = 0; i < numberOfXi; i++ )
	{
	flog( __FILE__, __LINE__ );
		switch( interpolationXi[i] )
		{
		case BASIS_LINEAR_LAGRANGE_INTERPOLATION:
	flog( __FILE__, __LINE__ );
			scaleFactorCount *= 2;
			err = fprintf( session->file, "l.Lagrange" );
			break;
		case BASIS_QUADRATIC_LAGRANGE_INTERPOLATION:
	flog( __FILE__, __LINE__ );
			scaleFactorCount *= 3;
             err = fprintf( session->file, "q.Lagrange" );
			break;
		case BASIS_CUBIC_LAGRANGE_INTERPOLATION:
	flog( __FILE__, __LINE__ );
			scaleFactorCount *= 4;
             err = fprintf( session->file, "c.Lagrange" );
			break;
		case BASIS_CUBIC_HERMITE_INTERPOLATION:
	flog( __FILE__, __LINE__ );
			scaleFactorCount *= 4;
             err = fprintf( session->file, "c.Hermite" );
			break;
		case BASIS_QUADRATIC1_HERMITE_INTERPOLATION:
	flog( __FILE__, __LINE__ );
			scaleFactorCount *= 4;
             err = fprintf( session->file, "q1.Hermite" );
			break;
		case BASIS_QUADRATIC2_HERMITE_INTERPOLATION:
	flog( __FILE__, __LINE__ );
			scaleFactorCount *= 4;
             err = fprintf( session->file, "q2.Hermite" );
			break;
		default:
	flog( __FILE__, __LINE__ );
			return FIELD_EXPORT_ERROR_UNKNOWN_INTERPOLATION;
		}

	flog( __FILE__, __LINE__ );
		if( ( err >= 0 ) && ( i < ( numberOfXi - 1 ) ) )
		{
			err = fprintf( session->file, "*" );
		}
		if( err < 0 )
		{
	flog( __FILE__, __LINE__ );
			return FIELD_EXPORT_ERROR_FILE_WRITE;
		}
	flog( __FILE__, __LINE__ );
	}

	flog( __FILE__, __LINE__ );
	switch( labelType )
	{
	case FIELD_IO_SCALE_FACTORS_NUMBER_TYPE:
		fprintf( session->file, ", #Scale factors=%d\n", scaleFactorCount );
		break;
	case FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE:
		fprintf( session->file, ", no modify, standard node based.\n" );
		break;
	default:
		return FIELD_EXPORT_ERROR_UNKNOWN_LABEL_TYPE;
	}
	flog( __FILE__, __LINE__ );

	return 0;
}


static char *FieldExport_GetVariableLabel( const int fieldType, const int variableType )
{
	switch( fieldType )
	{
	case FIELD_GEOMETRIC_TYPE:
		switch( variableType )
		{
		case FIELD_STANDARD_VARIABLE_TYPE:
			return "unknown";
		case FIELD_NORMAL_VARIABLE_TYPE:
            return "Normal_derivative,  field,  normal derivative of variable";
		case FIELD_TIME_DERIV1_VARIABLE_TYPE:
            return "first_time_derivative,  field,  first time derivative of variable";
		case FIELD_TIME_DERIV2_VARIABLE_TYPE:
            return "second_time_derivative,  field,  second time derivative of variable";
		default:
            return "unknown_geometry,  field,  unknown field variable type";
		}
	case FIELD_FIBRE_TYPE:
        switch( variableType )
		{
		case FIELD_STANDARD_VARIABLE_TYPE:
            return "fibres, anatomical, fibre";
		case FIELD_NORMAL_VARIABLE_TYPE:
            return "norm_der_fiber,  normal derivative of variable";
		case FIELD_TIME_DERIV1_VARIABLE_TYPE:
            return "first_time_fiber,  first time derivative of variable";
		case FIELD_TIME_DERIV2_VARIABLE_TYPE:
            return "second_time_fiber,  second time derivative of variable";
		default:
            return "unknown_fiber,  unknown field variable type";
		}
	case FIELD_GENERAL_TYPE:
        switch( variableType )
		{
		case FIELD_STANDARD_VARIABLE_TYPE:
			return "general,  field,  rectangular cartesian";
		case FIELD_NORMAL_VARIABLE_TYPE:
            return "norm_dev_variable,  field,  string";
		case FIELD_TIME_DERIV1_VARIABLE_TYPE:
            return "first_time_variable,  field,  first time derivative of variable";
		case FIELD_TIME_DERIV2_VARIABLE_TYPE:
            return "second_time_variable,  field,  second time derivative of variable";
		default:
            return "unknown_general,  field,  unknown field variable type";
		}
	case FIELD_MATERIAL_TYPE:
        switch( variableType )
		{
		case FIELD_STANDARD_VARIABLE_TYPE:
            return "material,  field,  rectangular cartesian";
		case FIELD_NORMAL_VARIABLE_TYPE:
            return "normal_material,  field,  normal derivative of variable";
		case FIELD_TIME_DERIV1_VARIABLE_TYPE:
            return "fist_time_material,  field,  first time derivative of variable";
		case FIELD_TIME_DERIV2_VARIABLE_TYPE:
            return "second_time_material,  field,  second time derivative of variable";
		default:
            return "unknown material,  field,  unknown field variable type";
		}
	default:
        switch( variableType )
		{
		case FIELD_STANDARD_VARIABLE_TYPE:
            return "unknown,  field,  unknown standand variable type";
		case FIELD_NORMAL_VARIABLE_TYPE:
            return "unknown,  field,  unknown normal derivative of variable";
		case FIELD_TIME_DERIV1_VARIABLE_TYPE:
            return "unknown,  field,  unknown first time derivative of variable";
		case FIELD_TIME_DERIV2_VARIABLE_TYPE:
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
	flog( __FILE__, __LINE__ );

	fblog( "FieldExport_GetCoordinateComponentLabel:coordinateSystem", coordinateSystem );

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

	flog( __FILE__, __LINE__ );
	return NULL;
}


static int FieldExport_File_NodeCount( const FileSession *const session, const int nodeCount )
{
	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	if( fprintf( session->file, " #Nodes=           %d\n", nodeCount ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_FieldCount( const FileSession *const session, const int fieldCount )
{
	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	if( fprintf( session->file, " #Fields=%d\n", fieldCount ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_CoordinateVariable( const FileSession *const session, const int variableIndex,
										CMISS_CoordinateSystem *coordinateSystem, const int componentCount )
{
	char *coordinateLabel;

	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	coordinateLabel = FieldExport_GetCoordinateVariableLabel( coordinateSystem );
	if( fprintf( session->file, " %d) %s, #Components=%d\n", variableIndex, coordinateLabel, componentCount ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_Variable( const FileSession *const session, const int variableIndex, const int fieldType, const int variableType, const int componentCount )
{
	char *variableLabel;

	if( session == NULL )
	{
		return FIELD_EXPORT_ERROR_BAD_HANDLE;
	}

	variableLabel = FieldExport_GetVariableLabel( fieldType, variableType );
	if( fprintf( session->file, " %d) %s, #Components=%d\n", variableIndex, variableLabel, componentCount ) < 0 )
	{
		return FIELD_EXPORT_ERROR_FILE_WRITE;
	}

	return 0;
}


static int FieldExport_File_CoordinateComponent( const FileSession *const session, CMISS_CoordinateSystem * coordinateSystem,
	const int componentNumber, const int numberOfXi, const int *const interpolationXi )
{
	const char * const componentLabel = FieldExport_GetCoordinateComponentLabel( coordinateSystem, componentNumber );

	flog( __FILE__, __LINE__ );

	fblog( "FieldExport_File_CoordinateComponent interpolationXi = ", interpolationXi );
	if( componentLabel == NULL )
	{
		fprintf( session->file, "   %d.  ", componentNumber );
	}
	else
	{
		fprintf( session->file, "   %s.  ", componentLabel );
	}

	flog( __FILE__, __LINE__ );

	return FieldExport_File_LagrangeHermiteScaleFactors( session, FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE, numberOfXi, interpolationXi );
}


static int FieldExport_File_Component( const FileSession *const session,
	const int componentNumber, const int numberOfXi, const int *const interpolationXi )
{
	flog( __FILE__, __LINE__ );

	fprintf( session->file, "   %d.  ", componentNumber );

	return FieldExport_File_LagrangeHermiteScaleFactors( session, FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE, numberOfXi, interpolationXi );
}


static int FieldExport_File_Nodes( const FileSession *const session, const int nodeCount, const int *const derivativeCount,
	const int *const elementDerivatives, int firstScaleIndex )
{
	int i, j;

	fprintf( session->file, "     #Nodes= %d\n", nodeCount );

	fblog( "derivatives = ", elementDerivatives );
	fblog( "derivatives[0] = ", elementDerivatives[0] );
	for( i = 0; i < nodeCount; i++ )
	{
		fprintf( session->file, " %5d.  #Values=%d\n", i + 1, elementDerivatives[i] );
		fprintf( session->file, "      Value indices:  " );
		for( j = 0; j < derivativeCount[i]; j++ )
		{
			fprintf( session->file, " %3d", elementDerivatives[ i + (j * nodeCount ) ] );
		}
		fprintf( session->file, "\n" );
		fprintf( session->file, "      Scale factor indices: " );
		for( j = 0; j < derivativeCount[i]; j++ )
		{
			firstScaleIndex++;
			fprintf( session->file, " %3d", firstScaleIndex );
		}
		fprintf( session->file, "\n" );
	}

	flog( __FILE__, __LINE__ );

	return 0;
}


static int FieldExport_File_ElementIndex( FileSession *session, const int dimensionCount, const int elementIndex )
{
	if( dimensionCount == 3 )
	{
          fprintf( session->file, " Element:            %d 0 0\n", elementIndex );
	}
	else if( dimensionCount == 2 )
	{
          fprintf( session->file, " Element:            0 %d 0\n", elementIndex );
	}
	else if( dimensionCount == 1 )
	{
          fprintf( session->file, " Element:            0 0 %d\n", elementIndex );
	}

	return 0;
}


static int FieldExport_File_ElementNodeIndices( FileSession *session, const int nodeCount, const int *const indices )
{
	int i;

	fprintf( session->file, "  Nodes:\n" );

	for( i = 0; i < nodeCount; i++ )
	{
		fprintf( session->file, "  %10d", indices[i] );
	}
	fprintf( session->file, "\n" );

	return 0;
}


static int FieldExport_File_ElementNodeScales( FileSession *session, const int scaleCount, const double *const scales )
{
	int i;

	fprintf( session->file, "   Scale factors:\n" );

	for( i = 0; i < scaleCount; i++ )
	{
		fprintf( session->file, "   %.16E", scales[i] );
	}
	fprintf( session->file, "\n" );

	return 0;
}


static FileSession *FieldExport_GetFileSession( const int handle )
{
	FileSessionListEntry *entry = fileSessions.next;

	while( entry != NULL )
	{
		if( entry->session->handle == handle )
		{
			return entry->session;
		}
		entry = entry->next;
	}

	return NULL;
}


static int FieldExport_OpenFileSession( const char *const name, int * const handle )
{
	int i;
	FileSessionListEntry *entry;
	FileSession *session = calloc( 1, sizeof( FileSession ) );

	session->handle = FieldExport_GetNewSessionHandle( EXPORT_TYPE_FILE );
	session->file = fopen( name, "w" );

	if( session->file == NULL )
	{
		free( session );
		return FIELD_EXPORT_ERROR_FILE_IO;
	}

	flog( "New session", session->handle );

	*handle = session->handle;

	entry = calloc( 1, sizeof( FileSessionListEntry ) );

	entry->session = session;
	entry->next = fileSessions.next;
	fileSessions.next = entry;

	return 0;
}


static int FieldExport_File_CloseSession( FileSession *session )
{
	FileSessionListEntry *entry = fileSessions.next;
	FileSessionListEntry *previousEntry = &fileSessions;

	handleTypes[session->handle] = EXPORT_TYPE_CLOSED;

	while( entry != NULL )
	{
		if( entry->session != session )
		{
			previousEntry = entry;
			entry = entry->next;
			continue;
		}

		fclose( entry->session->file );
		free( entry->session );

		previousEntry->next = entry->next;
		free( entry );
		break;
	}

	return 0;
}


static int FieldExport_File_NodeNumber( FileSession *session, const int nodeNumber )
{
	flog( __FILE__, __LINE__ );

	fprintf( session->file, " Node:            %d\n", nodeNumber );

	flog( __FILE__, __LINE__ );

	return 0;
}


static int FieldExport_File_NodeValues( FileSession *session, const int valueCount, const double *const values )
{
	int i;

	flog( __FILE__, __LINE__ );

	for( i = 0; i < valueCount; i++ )
	{
		fprintf( session->file, "  %.16E", values[i] );
	}
	fprintf( session->file, "\n" );

	flog( __FILE__, __LINE__ );

	return 0;
}


static void FieldExport_FieldDerivateLabels( FileSession *session, const int numberOfDerivatives, const int *const derivatives )
{
	int i;

	if( ( numberOfDerivatives == 1 ) && ( derivatives[0] == NO_PART_DERIV ) )
	{
		return;
	}

	for( i = 0; i < numberOfDerivatives; i++ )
	{
		switch( derivatives[i] )
		{
		case NO_PART_DERIV:
			break;
		case PART_DERIV_S1:
              fprintf( session->file, ", d/ds1" );
			  break;
		case PART_DERIV_S1_S1:
              fprintf( session->file, ", d2/ds1ds1" );
			  break;
		case PART_DERIV_S2:
              fprintf( session->file, ", d/ds2" );
			  break;
		case PART_DERIV_S2_S2:
              fprintf( session->file, ", d2/ds2ds2" );
			  break;
		case PART_DERIV_S1_S2:
              fprintf( session->file, ", d/ds3" );
			  break;
		case PART_DERIV_S3:
              fprintf( session->file, ", d2/ds3ds3" );
			  break;
		case PART_DERIV_S3_S3:
              fprintf( session->file, ", d2/ds3ds3" );
			  break;
		case PART_DERIV_S1_S3:
              fprintf( session->file, ", d2/ds1ds3" );
			  break;
		case PART_DERIV_S2_S3:
              fprintf( session->file, ", d2/ds2ds3" );
			  break;
		case PART_DERIV_S1_S2_S3:
              fprintf( session->file, ", d3/ds1ds2ds3" );
			  break;
		case PART_DERIV_S4:
              fprintf( session->file, ", d/ds4" );
			  break;
		case PART_DERIV_S4_S4:
              fprintf( session->file, ", d2/ds4ds4" );
			  break;
		case PART_DERIV_S1_S4:
              fprintf( session->file, ", d2/ds1ds4" );
			  break;
		case PART_DERIV_S2_S4:
              fprintf( session->file, ", d2/ds2ds4" );
			  break;
		case PART_DERIV_S3_S4:
              fprintf( session->file, ", d2/ds3ds4" );
			  break;
		case PART_DERIV_S1_S2_S4:
              fprintf( session->file, ", d3/ds1ds2ds4" );
			  break;
		case PART_DERIV_S1_S3_S4:
              fprintf( session->file, ", d3/ds1ds3ds4" );
			  break;
		case PART_DERIV_S2_S3_S4:
              fprintf( session->file, ", d3/ds2ds3ds4" );
			  break;
		case PART_DERIV_S1_S2_S3_S4:
              fprintf( session->file, ", d4/ds1ds2ds3ds4" );
			  break;
		default:
              fprintf( session->file, "unknown field variable type %d", derivatives[i] );
		}
	}
}


static int FieldExport_File_DerivativeIndices( FileSession *session, const int componentNumber, const int fieldType,
    const int variableType, const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
	flog( __FILE__, __LINE__ );

	//MUSTDO add a proper GetComponentLabel( fieldType, variableType, componentNumber ) function.
	fprintf( session->file, "   %d.  Value index= %d, #Derivatives= %d", componentNumber, valueIndex, numberOfDerivatives - 1 );

	FieldExport_FieldDerivateLabels( session, numberOfDerivatives, derivatives );

	fprintf( session->file, "\n" );

	flog( __FILE__, __LINE__ );

	return 0;
}


static int FieldExport_File_CoordinateDerivativeIndices( FileSession *session, const int componentNumber, CMISS_CoordinateSystem * coordinateSystem,
	const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
	const char * const componentLabel = FieldExport_GetCoordinateComponentLabel( coordinateSystem, componentNumber );

	flog( __FILE__, __LINE__ );

	//MUSTDO add a proper GetComponentLabel( fieldType, variableType, componentNumber ) function.
	if( componentLabel == NULL )
	{
		fprintf( session->file, "   %d.  Value index= %d, #Derivatives= %d", componentNumber, valueIndex, numberOfDerivatives - 1 );
	}
	else
	{
		fprintf( session->file, "   %s.  Value index= %d, #Derivatives= %d", componentLabel, valueIndex, numberOfDerivatives - 1 );
	}

	FieldExport_FieldDerivateLabels( session, numberOfDerivatives, derivatives );

	fprintf( session->file, "\n" );

	flog( __FILE__, __LINE__ );

	return 0;
}


/*
	Public API implementation
*/
int FieldExport_OpenSession( const int type, const char *const name, int * const handle )
{
	flog( __FILE__, __LINE__ );
	if( type == EXPORT_TYPE_FILE )
	{
		return FieldExport_OpenFileSession( name, handle );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_Group( const int handle, const char *const label )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}

	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_Group( fileSession, label );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_MeshDimensions( const int handle, const int dimensions )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_MeshDimensions( fileSession, dimensions );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_ScalingFactorCount( const int handle, const int scalingFactorCount )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_ScalingFactorCount( fileSession, scalingFactorCount );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_LagrangeHermiteScaleFactors( const int handle, const int numberOfXi, const int* const interpolationXi )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_LagrangeHermiteScaleFactors( fileSession, FIELD_IO_SCALE_FACTORS_NUMBER_TYPE, numberOfXi, interpolationXi );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_NodeCount( const int handle, const int nodeCount )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_NodeCount( fileSession, nodeCount );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_FieldCount( const int handle, const int fieldCount )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_FieldCount( fileSession, fieldCount );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_CoordinateVariable( const int handle, const int variableNumber, CMISS_CoordinateSystem coordinateSystem,
	const int componentCount )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_CoordinateVariable( fileSession, variableNumber, &coordinateSystem, componentCount );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_Variable( const int handle, const int variableNumber, const int fieldType, const int variableType,
	const int componentCount )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_Variable( fileSession, variableNumber, fieldType, variableType, componentCount );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_CoordinateComponent( const int handle, CMISS_CoordinateSystem coordinateSystem,
	const int componentNumber, const int numberOfXi, const int * const interpolationXi )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_CoordinateComponent( fileSession, &coordinateSystem, componentNumber, numberOfXi, interpolationXi );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_Component( const int handle, const int componentNumber, const int numberOfXi, const int * const interpolationXi )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_Component( fileSession, componentNumber, numberOfXi, interpolationXi );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_Nodes( const int handle, const int nodeCount, const int *const derivativeCount,
	const int *const elementDerivatives, const int firstScaleIndex )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_Nodes( fileSession, nodeCount, derivativeCount, elementDerivatives, firstScaleIndex );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_ElementIndex( const int handle, const int dimensionCount, const int index )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_ElementIndex( fileSession, dimensionCount, index );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_ElementNodeIndices( const int handle, const int nodeCount, const int* const indices )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_ElementNodeIndices( fileSession, nodeCount, indices );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_ElementNodeScales( const int handle, const int scaleCount, const double* const scales )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_ElementNodeScales( fileSession, scaleCount, scales );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}
}


int FieldExport_CloseSession( const int handle )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_CloseSession( fileSession );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}

	return 0;
}


int FieldExport_NodeNumber( const int handle, const int nodeNumber )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_NodeNumber( fileSession, nodeNumber );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}

	return 0;
}


int FieldExport_NodeValues( const int handle, const int valueCount, const double* const values )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_NodeValues( fileSession, valueCount, values );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}

	return 0;
}


int FieldExport_CoordinateDerivativeIndices( const int handle, const int componentNumber, CMISS_CoordinateSystem coordinateSystem,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_CoordinateDerivativeIndices( fileSession, componentNumber, &coordinateSystem, numberOfDerivatives, derivatives, valueIndex );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}

	return 0;
}


int FieldExport_DerivativeIndices( const int handle, const int componentNumber, const int fieldType, const int variableType,
    const int numberOfDerivatives, const int *const derivatives, const int valueIndex )
{
	int err = FieldExport_CheckHandle( handle );
	if( err != FIELD_EXPORT_NO_ERROR )
	{
		return err;
	}
	
	flog( __FILE__, __LINE__ );
	if( handleTypes[handle] == EXPORT_TYPE_FILE )
	{
		FileSession *fileSession = FieldExport_GetFileSession( handle );
		return FieldExport_File_DerivativeIndices( fileSession, componentNumber, fieldType, variableType, numberOfDerivatives, derivatives, valueIndex );
	}
	else
	{
		return FIELD_EXPORT_ERROR_UNKNOWN_TYPE;
	}

	return 0;
}
