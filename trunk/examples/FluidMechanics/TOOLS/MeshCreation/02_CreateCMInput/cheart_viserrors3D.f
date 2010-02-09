      SUBROUTINE PRINTZ(I,K)
!       Subprogram PRINTZ - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine prints messages and status for the preparatory routines
!
        USE CMGUI_VARS
        INTEGER I,K
        IF (I == 1) THEN
          PRINT *, ' ==>>  COMMAND LINE INPUTS COLLECTED  <<== '
        END IF
        IF (I == 2) PRINT *, ' ==>>  STATUS :: INPUT FILES READ  <<=='
        IF (I == 4) THEN
          PRINT *, ' ==>>  OUTPUT FILES  <<== '
        END IF
        IF (I == 100) PRINT *, ' ==>>  STATUS :: DONE  <<=='
      END SUBROUTINE PRINTZ


      SUBROUTINE ERRORZ(I)
!       Subprogram ERRORZ - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine prints error messages and stops execution for the preparatory routines
!
        USE CMGUI_VARS
        INTEGER I,K
        IF (I == 1) PRINT *, ' ==>>  FORCE QUIT :: INCOMPLETE ELEMENT TYPE (NO TRI, QUAD, TET, HEX SPECIFICATION)  <<== '
        IF (I == 2) PRINT *, ' ==>>  FORCE QUIT :: NO CALL LETTER or M BASIS  <<== '
        IF (I == 3) PRINT *, ' ==>>  FORCE QUIT :: NO FIELDS IN BASIS FILE <<== '
        IF (I == 4) PRINT *, ' ==>>  FORCE QUIT :: INCOMPLETE MAP or COEFFICIENT FILE(S) <<== '
        IF (I == 5) PRINT *, ' ==>>  FORCE QUIT :: PROBLEM ENCOUNTERED READING BASIS FILE <<== '
        IF (I == 6) PRINT *, ' ==>>  FORCE QUIT :: MULTIPLE ELEMENTS POINTING TO A SINGLE FACET <<== '
        PRINT *, ' ==>> QUITING <<== '
        STOP
      END SUBROUTINE ERRORZ


