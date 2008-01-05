/* \file
 * $Id: f90c_c.c 27 2007-07-24 16:52:51Z cpb $
 * \author Chris Bradley
 * \brief This file provides c utility routines for the f90_c module.
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
 * The Original Code is openCMISS
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

/*
File: f90_c.c
=================

This file provides c utility routines for the f90_c module.

$Id: f90c_c.c 27 2007-07-24 16:52:51Z cpb $

Functions included:

CStringLen           Returns the length of a c string
PackCharacters       Packs characters into a c string
UnpackCharacters     Unpacks characters from a c string

*/

/* Included files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Defines */

#ifdef VMS
#define CStringLen CSTRINGLEN
#define PackCharacters PACKCHARACTERS
#define UnPackCharacters UNPACKCHARACTERS
#endif
#ifdef mips
#define CStringLen cstringlen_
#define PackCharacters packcharacters_
#define UnPackCharacters unpackcharacters_
#endif
#ifdef WIN32
#define CStringLen cstringlen
#define PackCharacters packcharacters
#define UnPackCharacters unpackcharacters
#endif
#ifdef linux
#define CStringLen cstringlen
#define PackCharacters packcharacters
#define UnPackCharacters unpackcharacters
#endif
#ifdef _AIX
#define CStringLen cstringlen
#define PackCharacters packcharacters
#define UnPackCharacters unpackcharacters
#endif

/* Item Type defines */

/* Type definitions */

typedef int integer;
typedef integer logical;

/* Function prototypes */

void CStringLen(int *length,
  char *string);
void PackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string);
void UnPackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string);

/* Function Definitions */

/* Global variables */

/* Code */

void CStringLen(integer *length,
  char *string)
{
  *length=strlen(string);
}

void PackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string)
/*
Packs a fortran string (as a string of integers) into a c string
(as a string of characters).
*/

{
  char chr;

  chr = (char)*integer_char;
  *(integer_string+*char_num) = chr;

}

void UnPackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string)

/*
Unpacks a c string (as a string of characters) into a fortran string
(as a string of integers).
*/

{
  char chr;

  chr = *(integer_string+*char_num);
  *integer_char = (int)chr;

}
