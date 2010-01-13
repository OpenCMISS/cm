/* 
 * \file
 * $Id: opencmiss.h 582 2009-07-08 03:51:19Z catalept $
 * \author Chris Bradley
 * \brief The OpenCMISS library C header file.
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

#ifndef OPENCMISS_H
#define OPENCMISS_H
#endif

/* 
 * Defines
 */

const int CMISSNoError=0;

const int CMISSTrue=1;
const int CMISSFalse=0;

/*
 *=================================================================================================================================
 *
 * ANALYTIC_ANALYSIS_ROUTINES
 *
 *==================================================================================================================================
 */

/*
 *=================================================================================================================================
 *
 * BASE_ROUTINES
 *
 *==================================================================================================================================
 */

/* \addtogroup OPENCMISS_DiagnosticAndTimingConstants OPENCMISS::DiagnosticAndTiming::Constants
 * \brief Diagnostic and Timing constants.
 * @{  
 * \addtogroup OPENCMISS_DiagnosticTypes OPENCMISS::DiagnosticTypes
 * \brief Diganostic constants.
 * \see OPENCMISS::DiagnosticAndTiming,OPENCMISS
 * @}
 */  

#define CMISSAllDiagType 1
#define CMISSInDiagType 2
#define CMISSFromDiagType 3

/*
 *=================================================================================================================================
 *
 * CMISS
 *
 *==================================================================================================================================
 */

/* 
 * Struct defs
 */

struct CMISSCoordinateSystemType_;
struct CMISSRegionType_;

/* 
 * Type defs
 */



typedef int CMISSError;
typedef struct CMISSCoordinateSystemType_ *CMISSCoordinateSystemType;
typedef struct CMISSRegionType_ *CMISSRegionType;

/* 
 * Protypes
 */

/*
 *=================================================================================================================================
 *
 * Types
 *
 *==================================================================================================================================
 */

CMISSError CMISSCoordinateSystemTypeFinalise(CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSCoordinateSystemTypeInitialise(CMISSCoordinateSystemType *CoordinateSystem);

CMISSError CMISSRegionTypeFinalise(CMISSRegionType *Region);

CMISSError CMISSRegionTypeInitialise(CMISSRegionType *Err);

/*
 *=================================================================================================================================
 *
 * CMISS
 *
 *==================================================================================================================================
 */

CMISSError CMISSFinalise();

CMISSError CMISSInitialiseNum(int *WorldCoordinateSystemUserNumber,
			      int *WorldRegionUserNumber);

CMISSError CMISSInitialise(CMISSCoordinateSystemType *WorldCoordinateSystem,
			   CMISSRegionType *WorldRegion);

/*
 *=================================================================================================================================
 *
 * REGION
 *
 *==================================================================================================================================
 */

CMISSError CMISSRegionCreateFinish(CMISSRegionType Region);

CMISSError CMISSRegionCreateFinishNum(const int RegionUserNumber);

CMISSError CMISSRegionCreateStart(const int RegionUserNumber,
				  const CMISSRegionType ParentRegion,
				  CMISSRegionType Region);

CMISSError CMISSRegionCreateStartNum(const int RegionUserNumber,
				     const int ParentRegionUserNumber);

CMISSError CMISSRegionLabelGet(const CMISSRegionType Region,
			       const int LabelSize,
			       char *Label);

CMISSError CMISSRegionLabelGetNum(const int RegionUserNumber,
				  const int LabelSize,
				  char *Label);

CMISSError CMISSRegionLabelSet(const CMISSRegionType Region,
			       const int LabelSize,
			       const char *Label);

CMISSError CMISSRegionLabelSetNum(const int RegionUserNumber,
				  const int LabelSize,
				  const char *Label);
