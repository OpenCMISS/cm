!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all finite element base routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module contains all finite element base routines.
MODULE FINITE_ELEMENT_ROUTINES

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE STRINGS
  USE TYPES
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables
  
  !Interfaces

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE FEM_ELEMENT_MATRICES_FINALISE(ELEMENT_MATRICES,ERR,ERROR,*)

    !#### Subroutine: FEM_ELEMENT_MATRICES_FINALISE
    !###  Description:
    !###    Finalise the FEM element matrices for a problem and deallocate all memory.
    
    !Argument variables
    TYPE(FEM_ELEMENT_MATRICES_TYPE), POINTER :: ELEMENT_MATRICES
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FEM_ELEMENT_MATRICES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENT_MATRICES)) THEN
      
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FEM_ELEMENT_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("FEM_ELEMENT_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("FEM_ELEMENT_MATRICES_FINALISE")
    RETURN 1
  END SUBROUTINE FEM_ELEMENT_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE FEM_ELEMENT_MATRICES_INITIALISE(PROBLEM,ELEMENT_MATRICES,ERR,ERROR,*)

    !#### Subroutine: FEM_ELEMENT_MATRICES_INITIALISE 
    !###  Description:
    !###    Initialise the FEM element matrices for a problem.
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(FEM_ELEMENT_MATRICES_TYPE), POINTER :: ELEMENT_MATRICES
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FEM_ELEMENT_MATRICES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FEM_ELEMENT_MATRICES_INITIALISE")
    RETURN
999 CALL ERRORS("FEM_ELEMENT_MATRICES_INITIALISE",ERR,ERROR)
    CALL EXITS("FEM_ELEMENT_MATRICES_INITIALISE")
    RETURN 1
  END SUBROUTINE FEM_ELEMENT_MATRICES_INITIALISE
  

END MODULE FINITE_ELEMENT_ROUTINES
