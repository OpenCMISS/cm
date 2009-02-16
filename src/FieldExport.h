#ifndef FIELD_EXPORT_H
#define FIELD_EXPORT_H

int fieldexport_opensession( const int type, const char *const name, int * const handle );

int fieldexport_group( const int handle, const char *const label );

int fieldexport_meshdimensions( const int handle, const int dimensions );

int fieldexport_scalingfactorcount( const int handle, const int scalingFactorCount );

int fieldexport_lagrangehermitescalefactors( const int handle, const int numberOfXi, const int* const interpolationXi );

int fieldexport_nodecount( const int handle, const int nodeCount );

int fieldexport_fieldcount( const int handle, const int fieldCount );

int fieldexport_coordinatevariable( const int handle, const int variableNumber, CMISS_CoordinateSystem * coordinateSystem,
	const int componentCount );

int fieldexport_variable( const int handle, const int variableNumber, const int fieldType, const int variableType,
	const int componentCount );

int fieldexport_coordinatecomponent( const int handle, CMISS_CoordinateSystem * coordinateSystem,
	const int componentNumber, const int numberOfXi, const int * const interpolationXi );

int fieldexport_component( const int handle, const int componentNumber, const int numberOfXi, const int * const interpolationXi );

int fieldexport_nodes( const int handle, const int nodeCount, const int *const derivativeCount,
	const int *const elementDerivatives, const int firstScaleIndex );

int fieldexport_elementindex( const int handle, const int dimensionCount, const int index );

int fieldexport_elementnodeindices( const int handle, const int nodeCount, const int* const indices );

int fieldexport_elementnodescales( const int handle, const int scaleCount, const double* const scales );

int fieldexport_closesession( const int handle );

#endif