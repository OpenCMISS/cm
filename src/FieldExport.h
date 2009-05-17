/* \file
 * $Id: FieldExport.h 375 2009-02-27 11:05:03Z chrispbradley $
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
