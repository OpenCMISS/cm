/*
 * \file
 * $Id: MoreComplexMeshExample.c 20 2007-05-28 20:22:52Z cpb $
 * \author Chris Bradley
 * \brief This is an example program which sets up a field which uses a more complex mesh using OpenCMISS calls from C.
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
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss.h"

#define STRING_SIZE 20

#define REGION_USER_NUMBER 1

int main() 
{
  /* int WorldCoordinateSystemUserNumber;
     int WorldRegionUserNumber; */
  CMISSCoordinateSystemType *WorldCoordinateSystem=NULL;
  CMISSRegionType *WorldRegion=NULL,*Region=NULL;
  char Label[STRING_SIZE];
  int Err;

  /* Err = CMISSInitialiseNum(&WorldCoordinateSystemUserNumber,&WorldRegionUserNumber); */
  
  if(CMISSInitialise(&WorldCoordinateSystem,&WorldRegion) == CMISSNoError)
    {

      Err = CMISSRegionLabelGet(WorldRegion,STRING_SIZE,Label);
      printf("The world region label is '%s'.\n",Label);

      Err = CMISSRegionTypeInitialise(&Region);
      Err = CMISSRegionCreateStart(REGION_USER_NUMBER,WorldRegion,Region);
      Err = CMISSRegionLabelSet(Region,8,"Testing");
      Err = CMISSRegionCreateFinish(Region);

      Err = CMISSRegionLabelGet(Region,STRING_SIZE,Label);	       
      printf("The region label is '%s'.\n",Label);

      Err = CMISSRegionTypeFinalise(&Region);

      Err = CMISSFinalise();
    }

  return Err;
}
