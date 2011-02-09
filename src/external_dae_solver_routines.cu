/* \file
 * $Id: external_dae_solver_routines.c 1836 2010-12-20 17:25:14Z chrispbradley $
 * \author Chris Bradley
 * \brief This file provides the routines for solving differential-algebraic equations with an external solver.
 *.
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

/*
File: external_dae_solver_routines.c
===================
 
This file provides provides the routines for solving differential-algebraic equations with an external solver.

Functions included:

SolverDAEExternalIntegrate     Solves the differential-algebraic equation.

*/

/* Included files */

#include <stdio.h>
#include <stdlib.h>

#include "external_dae_solver_routines.h"
#include "cuda_solver_routines.cu"

/* Type definitions */

/* Function Definitions */
extern "C"
{
void SolverDAEExternalIntegrate(const int NumberOfDofs,
				const double StartTime,
				const double EndTime,
				double *InitialStep,
				const int OnlyOneModelIndex,
				int *ModelsData,
				int NumberOfState,
				double *StateData,
				int NumberOfParameters,
				double *ParametersData,
				int NumberOfIntermediate,
				double *IntermediateData,
				int *err)
{

	printf("start %f end %f steps %d\n", StartTime, EndTime, (int)((EndTime-StartTime)/InitialStep[0]));

	solve(StateData, StartTime, EndTime, (int)((EndTime-StartTime)/InitialStep[0]), NumberOfDofs, 1024, 10, 2, NULL);

}
}

