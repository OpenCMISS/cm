!> \file
!> $Id: cmiss_parmetis.f90 9 2007-05-15 13:52:02Z cpb $
!> \author Chris Bradley
!> \brief This module is a CMISS buffer module to the ParMETIS library.
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

!> This module is a CMISS buffer module to the ParMETIS library.
MODULE CMISS_PARMETIS
  
  USE BASE_ROUTINES
  USE KINDS
  USE ISO_VARYING_STRING
  
  IMPLICIT NONE
 
  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  INTERFACE

    SUBROUTINE ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ncon, nparts, tpwgts, ubvec, &
      & options, edgecut, part, comm)
#ifdef WIN32
      !DEC$ ATTRIBUTES C, reference, alias:'_ParMETIS_V3_PartKway' :: ParMETIS_V3_PartKway
#endif      
      USE KINDS
      INTEGER(INTG) :: vtxdist(*)
      INTEGER(INTG) :: xadj(*)
      INTEGER(INTG) :: adjncy(*)
      INTEGER(INTG) :: vwgt(*)
      INTEGER(INTG) :: adjwgt(*)
      INTEGER(INTG) :: wgtflag
      INTEGER(INTG) :: numflag
      INTEGER(INTG) :: ncon
      INTEGER(INTG) :: nparts
      REAL(SP) :: tpwgts(*)
      REAL(SP) :: ubvec(*)
      INTEGER(INTG) :: options(*)
      INTEGER(INTG) :: edgecut
      INTEGER(INTG) :: part(*)
      INTEGER(INTG) :: comm
    END SUBROUTINE PARMETIS_V3_PARTKWAY

    SUBROUTINE ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, &
      & ubvec, options, edgecut, part, comm)
#ifdef WIN32
      !DEC$ ATTRIBUTES C, reference, alias:'_ParMETIS_V3_PartMeshKway' :: ParMETIS_V3_PartMeshKway
#endif      
      USE KINDS
      INTEGER(INTG) :: elmdist(*)
      INTEGER(INTG) :: eptr(*)
      INTEGER(INTG) :: eind(*)
      INTEGER(INTG) :: elmwgt(*)
      INTEGER(INTG) :: wgtflag
      INTEGER(INTG) :: numflag
      INTEGER(INTG) :: ncon
      INTEGER(INTG) :: ncommonnodes
      INTEGER(INTG) :: nparts
      REAL(SP) :: tpwgts(*)
      REAL(SP) :: ubvec(*)
      INTEGER(INTG) :: options(*)
      INTEGER(INTG) :: edgecut
      INTEGER(INTG) :: part(*)
      INTEGER(INTG) :: comm      
    END SUBROUTINE ParMETIS_V3_PartMeshKway
    
  END INTERFACE

  PUBLIC PARMETIS_PARTKWAY,PARMETIS_PARTMESHKWAY
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Buffer routine to the ParMetis ParMETIS_V3_PartKway routine.
  SUBROUTINE PARMETIS_PARTKWAY(VERTEX_DISTANCE,XADJ,ADJNCY,VERTEX_WEIGHT,ADJ_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
    & NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR,ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: VERTEX_DISTANCE(:)
    INTEGER(INTG), INTENT(IN) :: XADJ(:)
    INTEGER(INTG), INTENT(IN) :: ADJNCY(:)
    INTEGER(INTG), INTENT(IN) :: VERTEX_WEIGHT(:)
    INTEGER(INTG), INTENT(IN) :: ADJ_WEIGHT(:)
    INTEGER(INTG), INTENT(IN) :: WEIGHT_FLAG
    INTEGER(INTG), INTENT(IN) :: NUM_FLAG
    INTEGER(INTG), INTENT(IN) :: NCON
    INTEGER(INTG), INTENT(IN) :: NUMBER_PARTS
    REAL(SP), INTENT(IN) :: TP_WEIGHTS(:)
    REAL(SP), INTENT(IN) :: UB_VEC(:)
    INTEGER(INTG), INTENT(IN) :: OPTIONS(:)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_EDGES_CUT
    INTEGER(INTG), INTENT(OUT) :: PARTITION(:)
    INTEGER(INTG), INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PARMETIS_PARTKWAY",ERR,ERROR,*999)

    CALL ParMETIS_V3_PartKway(VERTEX_DISTANCE,XADJ,ADJNCY,VERTEX_WEIGHT,ADJ_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
      & NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR)
    IF(ERR/=0) THEN
      CALL FLAG_ERROR("ParMetis error in ParMETIS_V3_PartKway",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PARMETIS_PARTKWAY")
    RETURN
999 CALL ERRORS("PARMETIS_PARTKWAY",ERR,ERROR)
    CALL EXITS("PARMETIS_PARTKWAY")
    RETURN 1
  END SUBROUTINE PARMETIS_PARTKWAY

  !
  !================================================================================================================================
  !

  !>Buffer routine to the ParMetis ParMETIS_V3_PartMeshKway routine.
  SUBROUTINE PARMETIS_PARTMESHKWAY(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDEX,ELEMENT_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
    & NUMBER_COMMON_NODES,NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR,ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: ELEMENT_DISTANCE(:)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_PTR(:)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_INDEX(:)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_WEIGHT(:)
    INTEGER(INTG), INTENT(IN) :: WEIGHT_FLAG
    INTEGER(INTG), INTENT(IN) :: NUM_FLAG
    INTEGER(INTG), INTENT(IN) :: NCON
    INTEGER(INTG), INTENT(IN) :: NUMBER_COMMON_NODES
    INTEGER(INTG), INTENT(IN) :: NUMBER_PARTS
    REAL(SP), INTENT(IN) :: TP_WEIGHTS(:)
    REAL(SP), INTENT(IN) :: UB_VEC(:)
    INTEGER(INTG), INTENT(IN) :: OPTIONS(:)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_EDGES_CUT
    INTEGER(INTG), INTENT(OUT) :: PARTITION(:)
    INTEGER(INTG), INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("PARMETIS_PARTMESHKWAY",ERR,ERROR,*999)
    
    CALL ParMETIS_V3_PartMeshKway(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDEX,ELEMENT_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
      & NUMBER_COMMON_NODES,NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR)
    
    IF(ERR/=0) THEN
      CALL FLAG_ERROR("ParMetis error in ParMETIS_V3_PartMeshKway",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PARMETIS_PARTMESHKWAY")
    RETURN
999 CALL ERRORS("PARMETIS_PARTMESHKWAY",ERR,ERROR)
    CALL EXITS("PARMETIS_PARTMESHKWAY")
    RETURN 1
  END SUBROUTINE PARMETIS_PARTMESHKWAY

  !
  !================================================================================================================================
  !
    
END MODULE CMISS_PARMETIS
