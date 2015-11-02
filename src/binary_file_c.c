/* \file
 * \author Chris Bradley
 * \brief This file provides c utility routines for the binary_file module
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
File: binary_file.c
===================
 
This file provides c utility routines for the binary_file module.

Functions included:

BinaryCloseFile      Close a binary file
BinaryOpenFile       Open a binary file
BinaryReadFile       Read data from a binary file
BinarySetFile        Sets the position of a binary file
BinarySkipFile       Skip bytes in a binary file
BinaryWriteFile      Write data to a binary file
IsBinaryFileOpen     Returns whether or not a binary file is open
IsEndBinaryFile      Returns whether or not at eof of a binary file

*/

/* Included files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef VMS
#define BinaryCloseFile BINARYCLOSEFILE
#define BinaryOpenFile BINARYOPENFILE
#define BinaryReadFile BINARYREADFILE
#define BinarySetFile BINARYSETFILE
#define BinarySkipFile BINARYSKIPFILE
#define BinaryWriteFile BINARYWRITEFILE
#define IsBinaryFileOpen ISBINARYFILEOPEN
#define IsEndBinaryFile ISENDBINARYFILE
#endif
#ifdef mips
#define BinaryCloseFile binaryclosefile_
#define BinaryOpenFile binaryopenfile_
#define BinaryReadFile binaryreadfile_
#define BinarySetFile binarysetfile_
#define BinarySkipFile binaryskipfile_
#define BinaryWriteFile binarywritefile_
#define IsBinaryFileOpen isbinaryfileopen_
#define IsEndBinaryFile isendbinaryfile_
#endif
#ifdef WIN32
#define BinaryCloseFile binaryclosefile
#define BinaryOpenFile binaryopenfile
#define BinaryReadFile binaryreadfile
#define BinarySetFile binarysetfile
#define BinarySkipFile binaryskipfile
#define BinaryWriteFile binarywritefile
#define IsBinaryFileOpen isbinaryfileopen
#define IsEndBinaryFile isendbinaryfile
#endif
#ifdef linux
#define BinaryCloseFile binaryclosefile
#define BinaryOpenFile binaryopenfile
#define BinaryReadFile binaryreadfile_
#define BinarySetFile binarysetfile
#define BinarySkipFile binaryskipfile
#define BinaryWriteFile binarywritefile_
#define IsBinaryFileOpen isbinaryfileopen
#define IsEndBinaryFile isendbinaryfile_
#endif

/* Item Type defines */
/*   These should be the same as those in constants.f90 */
#define INTEGERTYPE           1
#define SHORTINTTYPE          2
#define LONGINTTYPE           3
#define FLOATTYPE             4
#define DOUBLETYPE            5
#define QUADRUPLETYPE         6
#define CHARTYPE              7
#define LOGICALTYPE           8
#define COMPLEXTYPE           9
#define DOUBLECOMPLEXTYPE    10
#define QUADRUPLECOMPLEXTYPE 11

/* Binary file defines */
/*   These should be the same as those in binary_file.f90 */
#define MAXBINFILES 99
#define SAMEENDIAN 0
#define FLIPENDIAN 1

/* Type definitions */

typedef int logical;

/* Function prototypes */

void BinaryCloseFile(int *fileid,
  int *err, 
  char *error_string);
void BinaryOpenFile(int *fileid,
  char *filename,
  char* access_code,
  int *err,
  char *error_string);
void BinaryReadFile(int *fileid,
  int *endian,
  int *number_of_items, 
  int *item_type,
  char *data,
  int *err,
  char *error_string);
void BinarySetFile(int *fileid,
  int *set_code,
  int *err,
  char *error_string);
void BinarySkipFile(int *fileid,
  int *number_of_bytes, 
  int *err,
  char *error_string);
void BinaryWriteFile(int *fileid,
  int *endian,
  int *number_of_items, 
  int *item_type,
  char *data,
  int *err,
  char *error_string);
void IsBinaryFileOpen(int *fileid,
  int *returncode,
  int *err,
  char *error_string);
void IsEndBinaryFile(int *fileid,
  int *returncode,
  int *err,
  char *error_string);

/* Function Definitions */


/* Global variables */

FILE *binaryfiles[MAXBINFILES]=
  {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

/* Code */

void BinaryCloseFile(int *fileid,
  int *err,
  char *error_string)

/*
Closes the binary file specified by fileid.
*/

{
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    if(!binaryfiles[*fileid-1])
    {
      *err = 0;      
    }
    else
    {
      *err = fclose(binaryfiles[*fileid-1]);
      binaryfiles[*fileid-1] = (FILE *)NULL;
      if(*err != 0)
      {
        strcpy(error_string,">>ERROR: error closing binary file");
      }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinaryOpenFile(int *fileid,
  char *filename,
  char *access_code,
  int *err,
  char *error_string)

/* 
Opens a binary file specified by fileid and name filename.
*/

{
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    if(binaryfiles[*fileid-1])
    {
      *err=1;
      strcpy(error_string,">>ERROR: binary file is already open");
    }
    else
    {
      if(binaryfiles[*fileid-1] = fopen(filename,access_code))
	    {
	      *err=0;
	    }
      else
	    {
	      *err=1;
	      strcpy(error_string,">>ERROR: binary file could not be opened");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinaryReadFile(int *fileid,
  int *endian,
  int *number_of_items, 
  int *item_type,
  char *data,
  int *err,
  char *error_string)

/* 
Reads number_of_items of data of a type given by item_type from
a binary file specified by fileid into an array iven by data.

The default endian ordering is big endian. This is specified by
endian=0. If little endian is required endian must be set to a number
other than 0.

*/

{
  int i,item_size,j,number_of_bytes,start_byte,temp;
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      if(INTEGERTYPE == *item_type)
	    {
	      /* Integer data */
	      item_size=sizeof(int);
	    }
      else if(FLOATTYPE == *item_type)
	    {
	      /* Float data */
	      item_size=sizeof(float);
	    }
      else if(DOUBLETYPE == *item_type)
	    /* Double data */
	    {
	      item_size=sizeof(double);
	    }
      else if(CHARTYPE == *item_type)
	    /* Character data */
	    {
	      item_size=sizeof(char);
	    }
      else if(LOGICALTYPE == *item_type)
	    /* Logical data */
	    {
	      item_size=sizeof(logical);
	    }
      else if(SHORTINTTYPE == *item_type)
	    /* Short integer data */
	    {
	      item_size=sizeof(short int);
	    }
      else
	    {
	      *err=1;
	      strcpy(error_string,">>ERROR: Invalid item type");
	    }
	  
      if(SAMEENDIAN == *endian || CHARTYPE == *item_type) 
	    {
	      /* Default big endian format */
	      number_of_bytes=*number_of_items * item_size;
	      for(i = 0; i < number_of_bytes ; i++)
        {
          temp = getc(binaryfile); 
          data[i] = (char)temp;  
        }
        if(CHARTYPE == *item_type)
        {
          data[number_of_bytes]='\0';
        }
	    }
      else
	    {
	      /* Little endian format - must reverse byte ordering */
	      for(i = 0; i < *number_of_items; i++)
        {
          start_byte=i*item_size;
          for(j = 0; j < item_size; j++)
          {
            temp = getc(binaryfile); 
            data[start_byte+item_size-j-1] = (char)temp;
          }
        }
	    }
	  
      *err=ferror(binaryfile);
      if(*err != 0)
	    {
	      strcpy(error_string,">>ERROR: error reading binary file");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinarySetFile(int *fileid,
  int *set_code,
  int *err,
  char *error_string)

/* 
Sets the position of the file pointer of the binary file (given
by fileid) to either the beginning (set_code=0), current
position (set_code=1) or end (set_code=2) of the file.
*/

{
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      switch(*set_code)
      {
        case 0: /* Beginning of a file */
        {
          *err=fseek(binaryfile,(long)0,SEEK_SET);
          if(*err != 0)
          {
            strcpy(error_string,">>ERROR: Could not set to beginning of file");
          }
        } break;
        case 1: /* Current file position */
        {
          *err=fseek(binaryfile,(long)0,SEEK_CUR);
          if(*err != 0)
          {
            strcpy(error_string,">>ERROR: Could not set to current file position");
          }
        } break;
        case 2: /* End of file */
        {
          *err=fseek(binaryfile,(long)0,SEEK_END);
          if(*err != 0)
          {
            strcpy(error_string,">>ERROR: Could not set to end of file position");
          }
        } break;
        default:
        {
          *err=1;
          strcpy(error_string,">>ERROR: Invalid set_code");          
        } break;
      }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinarySkipFile(int *fileid,
  int *number_of_bytes, 
  int *err,
  char *error_string)

/* 
Skips number_of_bytes bytes of data in a binary file specified
by fileid.
*/

{
  int i,temp;
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      for(i = 0; i < *number_of_bytes ; i++)
	    {
	      temp = getc(binaryfile); 
	    }
	  
      *err=ferror(binaryfile);
      if(*err != 0)
	    {
	      strcpy(error_string,">>ERROR: error skipping binary file");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinaryWriteFile(int *fileid,
  int *endian,
  int *number_of_items, 
  int *item_type,
  char *data,
  int *err,
  char *error_string)

/* 
Writes number_of_items of data of a type given byitem_type to a
binary file specified by fileid from an array given by data. 

The default endian ordering is big endian. This is specified by
endian=0. If little endian is required endian must be set to a number
other than 0.

*/

{
  int i,item_size,j,number_of_bytes,start_byte;
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      if(INTEGERTYPE == *item_type)
	    {
	      /* Integer data */
	      item_size=sizeof(int);
	    }
      else if(FLOATTYPE == *item_type)
	    {
	      /* Float data */
	      item_size=sizeof(float);
	    }
      else if(DOUBLETYPE == *item_type)
	    /* Double data */
	    {
	      item_size=sizeof(double);
	    }
      else if(CHARTYPE == *item_type)
	    /* Character data */
	    {
	      item_size=sizeof(char);
	    }
      else if(LOGICALTYPE == *item_type)
	    /* Logical data */
	    {
	      item_size=sizeof(logical);
	    }
      else if(SHORTINTTYPE == *item_type)
	    /* Short integer data */
	    {
	      item_size=sizeof(short int);
	    }
      else
	    {
	      *err=1;
	      strcpy(error_string,">>ERROR: Invalid item type");
	    }
	  
      if(SAMEENDIAN == *endian || CHARTYPE == *item_type) 
	    {
	      /* Default big endian format */
	      number_of_bytes=*number_of_items * item_size;
	      for(i = 0; i < number_of_bytes ; i++)
        {
          putc(data[i],binaryfile);
        }
	    }
      else
	    {
	      /* Little endian format - must reverse byte ordering */
	      for(i = 0; i < *number_of_items; i++)
        {
          start_byte=i*item_size;
          for(j = 0; j < item_size; j++)
          {
            putc(data[start_byte+item_size-j-1],binaryfile);
          }
        }
	    }
      
      *err=ferror(binaryfile);
      if(*err != 0)
	    {
	      strcpy(error_string,">>ERROR: error writing binary file");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void IsBinaryFileOpen(int *fileid,
  int *returncode,
  int *err,
  char *error_string)

/* 
Returns 1 in returncode if a binary file specified  by fileid
is open, 0 if not.
*/

{
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {      
      *returncode = 0;
    }
    else
    {
      *returncode = 1;
    }
    *err=0;
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}


void IsEndBinaryFile(int *fileid,
  int *returncode,
  int *err,
  char *error_string)

/* 
Returns 1 in returncode if a binary file specified by fileid
is at end of file (eof), 0 if not.
*/

{
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(feof(binaryfile))
    {      
      *returncode = 1;
    }
    else
    {
      *returncode = 0;
    }
    *err=0;
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

