/* \file
 * $Id: cmiss.c 1669 2010-10-27 17:54:14Z chrispbradley $
 * \author Chris Bradley
 * \brief This file contains system level routines.
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

/* Included files */

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

/* Type definitions */

/* Function prototypes */

void CMISSResetFatalHandler(void);

void CMISSSetFatalHandler(void);

void CMISSInitFatalHandler(void);

/* Internal functions */

static void CMISSFatalHandler(int sig,
#  if defined (sun)
                         siginfo_t *sip,
                         ucontext_t *uap);
#  else
			 int code,
			 struct sigcontext *sc);
#  endif

/* Static variables */

/* static sigjmp_buf jump_buffer; */ 
static struct sigaction fatal_sigaction;
static struct sigaction old_SIGBUS_action;
#ifdef SIGEMT
static struct sigaction old_SIGEMT_action;
#endif
static struct sigaction old_SIGFPE_action;
static struct sigaction old_SIGILL_action;
static struct sigaction old_SIGINT_action;
static struct sigaction old_SIGABRT_action;
static struct sigaction old_SIGSEGV_action;
static struct sigaction old_SIGTRAP_action;

void CMISSResetFatalHandler()
{
#if defined (SIGBUS)
  if( 0 != sigaction(SIGBUS,&old_SIGBUS_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGBUS handler.");
    }
#endif /* defined (SIGBUS) */
#ifdef SIGEMT
  if( 0 != sigaction(SIGEMT,&old_SIGEMT_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGEMT handler.");
    }
#endif
  if( 0 != sigaction(SIGFPE,&old_SIGFPE_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGFPE handler.");
    }
  if( 0 != sigaction(SIGILL,&old_SIGILL_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGILL handler.");
    }
  if( 0 != sigaction(SIGINT,&old_SIGINT_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGINT handler.");
    }
  if( 0 != sigaction(SIGABRT,&old_SIGABRT_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGABRT handler.");
    }
  if( 0 != sigaction(SIGSEGV,&old_SIGSEGV_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGSEGV handler.");
    }
#if defined (SIGTRAP)
  if( 0 != sigaction(SIGTRAP,&old_SIGTRAP_action,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not reset SIGTRAP handler.");
    }
#endif /* defined (SIGTRAP) */
}

void CMISSSetFatalHandler(void)
{
#if defined (unix) || defined (_AIX)
#if defined (SIGBUS)
  if( 0 != sigaction(SIGBUS,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: Could not set SIGBUS handler.");
    }
#endif /* defined (SIGBUS) */
#ifdef SIGEMT
  if( 0 != sigaction(SIGEMT,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGEMT handler.");
    }
#endif
  if( 0 != sigaction(SIGFPE,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGFPE handler.");
    }
  if( 0 != sigaction(SIGILL,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGILL handler.");
    }
  if( 0 != sigaction(SIGINT,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGINT handler.");
    }
  if( 0 != sigaction(SIGABRT,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGABRT handler.");
    }
  if( 0 != sigaction(SIGSEGV,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGSEGV handler.");
    }
#if defined (SIGTRAP)
  if( 0 != sigaction(SIGTRAP,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,">>WARNING: could not set SIGTRAP handler.");
    }
#endif /* defined (SIGTRAP) */
#endif /* defined (unix) || defined (_AIX) */
}

static void CMISSFatalHandler(int sig,
#  if defined (sun)
                         siginfo_t *sip,
                         ucontext_t *uap)
#  else
			 int code,
			 struct sigcontext *sc)
#  endif
{

#if defined(_AIX)
  /* this from libxlf90.a provides a good description of what went wrong */
  xl__sigdump(sig,code,sc);
#else
  switch(sig)
    {
#if defined (SIGBUS)
    case SIGBUS:
      {
	fprintf(stderr,">>FATAL ERROR: Bus error occurred.\n");
     } break;
#endif /* defined (SIGBUS) */
#if defined (SIGEMT)
    case SIGEMT:
      {
	fprintf(stderr,">>FATAL ERROR: EMT occurred.\n");
      } break;
#endif /* defined (SIGEMT) */
    case SIGFPE:
      {
	fprintf(stderr,">>FATAL ERROR: Floating point execption occurred.\n");
      } break;
    case SIGILL:
      {
	fprintf(stderr,">>FATAL ERROR: Illegal instruction occurred.\n");
      } break;
    case SIGINT:
      {
	fprintf(stderr,">>FATAL ERROR: Interrupt occurred.\n");
      } break;
    case SIGABRT:
      {
	fprintf(stderr,">>FATAL ERROR: Abort occurred.\n");
      } break;
    case SIGSEGV:
      {
	fprintf(stderr,">>FATAL ERROR: Segment violation occurred.\n");
      } break;
#if defined (SIGTRAP)
    case SIGTRAP:
      {
	fprintf(stderr,">>FATAL ERROR: Trace trap occurred.\n");
      } break;
#endif /* defined (SIGTRAP) */
    default:
      {
	fprintf(stderr,">>FATAL ERROR: Unknown signal %d occurred.\n",code);
      } break;
    }
#endif

  /* There is an issue with signal handling in a library such as OpenCMISS. The best option would be to 
     jump back to where the user called an OpenCMISS routine and let them process the error if that is 
     what they wish. This, however, would require setting the long jump buffer at each entry point to the 
     OpenCMISS library. This could be done by modifying enters in opencmiss.f90 but may cause performance
     problems (probably not too bad as the major computations are inside the library rather than at the 
     interface). For now just stop the program on a signal. */

  /* siglongjmp(jump_buffer,sig); */

  exit(0);
}

void CMISSInitFatalHandler(void)
{
  fatal_sigaction.sa_flags = SA_NODEFER;
  fatal_sigaction.sa_handler = (void (*)(int))CMISSFatalHandler;
  if( 0 != sigemptyset(&fatal_sigaction.sa_mask) )
    {
      fprintf(stderr,">>WARNING: sigemptyset failed in CMISSInitFatalHandler.");
    }

#if defined (SIGBUS)
  sigaction(SIGBUS,NULL,&old_SIGBUS_action);
#endif /* defined (SIGBUS) */
#if defined (SIGEMT)
  sigaction(SIGEMT,NULL,&old_SIGEMT_action);
#endif /* defined (SIGEMT) */
  sigaction(SIGFPE,NULL,&old_SIGFPE_action);
  sigaction(SIGILL,NULL,&old_SIGILL_action);
  sigaction(SIGINT,NULL,&old_SIGINT_action);
  sigaction(SIGABRT,NULL,&old_SIGABRT_action);
  sigaction(SIGSEGV,NULL,&old_SIGSEGV_action);
#if defined (SIGTRAP)
  sigaction(SIGTRAP,NULL,&old_SIGTRAP_action);
#endif /* defined (SIGTRAP) */
}

