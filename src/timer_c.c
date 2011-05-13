/* \file
 * $Id$
 * \author Chris Bradley
 * \brief This file provides c utility routines for the timer module.
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
File: timer_c.c
===================
 
This file provides c utility routines for the timer module.

Functions included:

CPUTimer      Returns the CPU time

*/

/* Included files */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef VMS
#include <time.h>
#include <types.h>
#endif
#ifdef mips
#include <sys/times.h>
#include <sys/types.h>
#endif
#ifdef WIN32
#include <time.h>
#endif
#ifdef linux
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#endif
#ifdef _AIX
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#endif

/* Type definitions */

/* Function Definitions */

void CPUTimer(double *return_time,
  int *flag,
  int *err,
  char error_string[])

/* keep these up to date with timer_f.f90 */

#define USER_CPU   1
#define SYSTEM_CPU 2
#define TOTAL_CPU  3

{
  double system_time,total_time,user_time;
#if defined VMS
  struct tbuffer current_time;

  times(&current_time);
  user_time = ((double)(current_time.proc_user_time))/100.0;
  system_time = ((double)(current_time.proc_system_time))/100.0;
#endif
#if defined mips
  struct tms current_time;

  times(&current_time);
  user_time = ((double)(current_time.tms_utime))/((double)(CLOCKS_PER_SEC));
  system_time = ((double)(current_time.tms_stime))/((double)(CLOCKS_PER_SEC));
#endif
#if defined WIN32
  clock_t current_time;

  current_time = clock();
  user_time = ((double)(current_time)/(double)(CLOCKS_PER_SEC));
  system_time = 0.0;
#endif
#if defined linux
  struct tms current_time;

  times(&current_time);
  user_time = ((double)(current_time.tms_utime))/((double)sysconf(_SC_CLK_TCK));
  system_time = ((double)(current_time.tms_stime))/((double)sysconf(_SC_CLK_TCK));
#endif
#if defined _AIX
  struct tms current_time;

  times(&current_time);
  user_time = ((double)(current_time.tms_utime))/((double)sysconf(_SC_CLK_TCK));
  system_time = ((double)(current_time.tms_stime))/((double)sysconf(_SC_CLK_TCK));
#endif

  total_time = user_time+system_time;
  *err=0;
  switch(*flag)
  {
    case TOTAL_CPU:
    {
      *return_time = total_time;
    }; break;
    case USER_CPU:
    {
      *return_time = user_time;
    }; break;
    case SYSTEM_CPU:
    {
      *return_time = system_time;
    }; break;
    default:
    {
      *err=1;
      strcpy(error_string,"Invalid operation code");
      *return_time = -99999.0; /* get some attention! */
    }
  }     
}


