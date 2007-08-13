!> \file
!> $Id: generated_mesh_routines.f90 9 2007-05-15 13:52:02Z cpb $
!> \author Chris Bradley
!> \brief This module handles all generated mesh routines.
!> \todo Move generated regular mesh from mesh routines and generalise.
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

!> This module handles all generated mesh routines.
MODULE GENERATED_MESH_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup GENERATED_MESH_ROUTINES_GeneratedMeshTypes GENERATED_MESH_ROUTINES::GeneratedMeshTypes 
  !> \brief Generated mesh types.
  !> \see GENERATED_MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_TYPE=1 !<A regular generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_POLAR_TYPE=2 !<A polar generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_FRACTAL_TREE_TYPE=3 !<A fractal tree generated mesh. \see GENERATED_MESH_ROUTINES_GeneratedMeshTypes,GENERATED_MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_CREATE_FINISH
    !###  Description:
    !###    Finishes the creation of a generated mesh.

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE GENERATED_MESH_CREATE_START(USER_NUMBER,REGION,GENERATED_MESH,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_CREATE_START
    !###  Description:
    !###    Starts the creation of a generated mesh.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("GENERATED_MESH_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(GENERATED_MESH)) THEN
        CALL FLAG_ERROR("Generated mesh is already associated",ERROR,*999)
      ELSE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_START")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_START",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_START")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE GENERATED_MESH_DESTROY(GENERATED_MESH,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_DESTROY
    !###  Description:
    !###    Destroys a generated mesh.

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE GENERATED_MESH_FINALISE(GENERATED_MESH,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_FINALISE
    !###  Description:
    !###    Finalises a generated mesh and dellocates all memory.

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("GENERATED_MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
    ENDIF
 
    CALL EXITS("GENERATED_MESH_FINALISE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_FINALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_FINALISE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_FINALISE

  !
  !================================================================================================================================
  !
  SUBROUTINE GENERATED_MESH_INITALISE(GENERATED_MESH,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_INITALISE
    !###  Description:
    !###    Initialises a generated mesh.

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("GENERATED_MESH_INITALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
      GENERATED_MESH%USER_NUMBER=0
      NULLIFY(GENERATED_MESH%REGION)
      GENERATED_MESH%TYPE=0
      NULLIFY(GENERATED_MESH%REGULAR_MESH)
      NULLIFY(GENERATED_MESH%MESH)
    ELSE
      CALL FLAG_ERROR("Generated mesh is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_INITALISE")
    RETURN
999 CALL ERRORS("GENERATED_MESH_INITALISE",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_INITALISE")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_INITALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE GENERATED_MESH_TYPE_SET(GENERATED_MESH,GENERATED_MESH_TYPE,ERR,ERROR,*)

    !#### Subroutine: GENERATED_MESH_TYPE_SET
    !###  Description:
    !###    Sets/changes the type of a generated mesh.

    !Argument variables
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    INTEGER(INTG), INTENT(IN) :: GENERATED_MESH_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(GENERATED_MESH)) THEN
    ELSE
      CALL FLAG_ERROR("Generated mesh is already associated",ERROR,*999)
    ENDIF
 
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN
999 CALL ERRORS("GENERATED_MESH_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GENERATED_MESH_CREATE_FINISH")
    RETURN 1   
  END SUBROUTINE GENERATED_MESH_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_CREATE_REGULAR(USER_NUMBER,REGION,ORIGIN,MAXIMUM_EXTENT,NUMBER_ELEMENTS_XI,BASIS,MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_CREATE_REGULAR
    !###  Description:
    !###    Creates the regular mesh with the given USER_NUMBER in the specifed REGION. The mesh starts at the ORIGIN(:) and has
    !###    a maximum extent position of MAXIMUM_EXTENT(:) with the NUMBER_OF_ELEMENTS(:) in each direction. Each element is of
    !###    the specified BASIS type. A pointer to the finished mesh is returned in MESH.  

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    REAL(DP), INTENT(IN) :: ORIGIN(:),MAXIMUM_EXTENT(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_ELEMENTS_XI(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: ni,ne,ne1,ne2,ne3,NN,nn1,nn2,nn3,np,np1,np2,np3,TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_ELEMENTS_XI(3), &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    REAL(DP) :: INITIAL_POSITION(3),DELTA_COORD(3),MY_ORIGIN(3),MY_EXTENT(3),MESH_SIZE(3)
    LOGICAL :: BASIS_OK,ELEMENTS_OK
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(NODES_TYPE), POINTER :: NODES

    CALL ENTERS("MESH_CREATE_REGULAR",ERR,ERROR,*999)

    NULLIFY(MESH)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
        !Determine the coordinate system and create the regular mesh for that system
        SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          NUMBER_OF_DIMENSIONS=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
          IF(SIZE(ORIGIN)==NUMBER_OF_DIMENSIONS) THEN
            IF(SIZE(MAXIMUM_EXTENT)==NUMBER_OF_DIMENSIONS) THEN
              IF(ASSOCIATED(BASIS)) THEN
                SELECT CASE(BASIS%TYPE)
                CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                  IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                    ELEMENTS_OK=ALL(NUMBER_ELEMENTS_XI>0)
                    IF(.NOT.ELEMENTS_OK) CALL FLAG_ERROR("Must have 1 or more elements in all directions",ERR,ERROR,*999)
                    BASIS_OK=ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED)
                    IF(.NOT.BASIS_OK) CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented",ERR,ERROR,*999)
                    !Calculate sizes
                    TOTAL_NUMBER_OF_NODES=1
                    TOTAL_NUMBER_OF_ELEMENTS=1
                    TOTAL_NUMBER_OF_NODES_XI=1
                    TOTAL_NUMBER_ELEMENTS_XI=0
                    DELTA_COORD=0.0_DP
                    DO ni=1,BASIS%NUMBER_OF_XI
                      TOTAL_NUMBER_OF_NODES_XI(ni)=(BASIS%NUMBER_OF_NODES_XI(ni)-2)*NUMBER_ELEMENTS_XI(ni)+ &
                        & NUMBER_ELEMENTS_XI(ni)+1
                      TOTAL_NUMBER_ELEMENTS_XI(ni)=NUMBER_ELEMENTS_XI(ni)
                      TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES*TOTAL_NUMBER_OF_NODES_XI(ni)
                      TOTAL_NUMBER_OF_ELEMENTS=TOTAL_NUMBER_OF_ELEMENTS*TOTAL_NUMBER_ELEMENTS_XI(ni)
                    ENDDO !ni
                    !Create the default node set
                    CALL NODES_CREATE_START(TOTAL_NUMBER_OF_NODES,REGION,NODES,ERR,ERROR,*999)
                    MY_ORIGIN=0.0_DP
                    MY_EXTENT=0.0_DP
                    MY_ORIGIN(1:NUMBER_OF_DIMENSIONS)=ORIGIN
                    MY_EXTENT(1:NUMBER_OF_DIMENSIONS)=MAXIMUM_EXTENT
                    MESH_SIZE=MY_EXTENT-MY_ORIGIN
                    DO ni=1,BASIS%NUMBER_OF_XI
                      !This assumes that the xi directions are aligned with the coordinate directions
                      DELTA_COORD(ni)=MESH_SIZE(ni)/REAL(NUMBER_ELEMENTS_XI(ni),DP)
                    ENDDO !ni
                    DO np3=1,TOTAL_NUMBER_OF_NODES_XI(3)
                      DO np2=1,TOTAL_NUMBER_OF_NODES_XI(2)
                        DO np1=1,TOTAL_NUMBER_OF_NODES_XI(1)
                          np=np1+(np2-1)*TOTAL_NUMBER_OF_NODES_XI(1)+(np3-1)*TOTAL_NUMBER_OF_NODES_XI(2)
                          INITIAL_POSITION(1)=MY_ORIGIN(1)+REAL(np1-1,DP)*DELTA_COORD(1)
                          INITIAL_POSITION(2)=MY_ORIGIN(2)+REAL(np2-1,DP)*DELTA_COORD(2)
                          INITIAL_POSITION(3)=MY_ORIGIN(3)+REAL(np3-1,DP)*DELTA_COORD(3)
                          CALL NODE_INITIAL_POSITION_SET(np,INITIAL_POSITION(1:NUMBER_OF_DIMENSIONS),NODES,ERR,ERROR,*999)
                        ENDDO !np1
                      ENDDO !np2
                    ENDDO !np3
                    CALL NODES_CREATE_FINISH(REGION,ERR,ERROR,*999)
                    !Create the mesh
                    CALL MESH_CREATE_START(USER_NUMBER,REGION,SIZE(NUMBER_ELEMENTS_XI,1),MESH,ERR,ERROR,*999)
                    !Create the elements
                    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
                    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,1,BASIS,MESH_ELEMENTS,ERR,ERROR,*999)
                    !Set the elements for the regular mesh
                    ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element nodes",ERR,ERROR,*999)                
                    !Step in the xi(3) direction
                    DO ne3=1,TOTAL_NUMBER_ELEMENTS_XI(3)+1
                      DO ne2=1,TOTAL_NUMBER_ELEMENTS_XI(2)+1
                        DO ne1=1,TOTAL_NUMBER_ELEMENTS_XI(1)+1
                          IF(BASIS%NUMBER_OF_XI<3.OR.ne3<=TOTAL_NUMBER_ELEMENTS_XI(3)) THEN
                            IF(BASIS%NUMBER_OF_XI<2.OR.ne2<=TOTAL_NUMBER_ELEMENTS_XI(2)) THEN
                              IF(ne1<=TOTAL_NUMBER_ELEMENTS_XI(1)) THEN
                                ne=ne1
                                np=1+(ne1-1)*(BASIS%NUMBER_OF_NODES_XI(1)-1)
                                IF(BASIS%NUMBER_OF_XI>1) THEN
                                  ne=ne+(ne2-1)*TOTAL_NUMBER_ELEMENTS_XI(1)
                                  np=np+(ne2-1)*TOTAL_NUMBER_OF_NODES_XI(1)*(BASIS%NUMBER_OF_NODES_XI(2)-1)
                                  IF(BASIS%NUMBER_OF_XI>2) THEN
                                    ne=ne+(ne3-1)*TOTAL_NUMBER_ELEMENTS_XI(1)*TOTAl_NUMBER_ELEMENTS_XI(2)
                                    np=np+(ne3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)* &
                                      & (BASIS%NUMBER_OF_NODES_XI(3)-1)
                                  ENDIF
                                ENDIF
                                nn=0
                                DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                  nn=nn+1
                                  ELEMENT_NODES(nn)=np+(nn1-1)                              
                                ENDDO !nn1
                                IF(BASIS%NUMBER_OF_XI>1) THEN
                                  DO nn2=2,BASIS%NUMBER_OF_NODES_XI(2)
                                    DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                      nn=nn+1
                                      ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)
                                    ENDDO !nn1
                                  ENDDO !nn2
                                  IF(BASIS%NUMBER_OF_XI>2) THEN
                                    DO nn3=2,BASIS%NUMBER_OF_NODES_XI(3)
                                      DO nn2=1,BASIS%NUMBER_OF_NODES_XI(2)
                                        DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                          nn=nn+1
                                          ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)+ &
                                            & (nn3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                        ENDDO !nn1
                                      ENDDO !nn2
                                    ENDDO !nn3
                                  ENDIF
                                ENDIF
                                CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS,ELEMENT_NODES, &
                                  & ERR,ERROR,*999)
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO !ne1
                      ENDDO !ne2
                    ENDDO !ne3
                    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,1,ERR,ERROR,*999)
                    !Finish the mesh
                    CALL MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("The number of xi directions of the given basis does not match the size of &
                      &the number of elements for the mesh",ERR,ERROR,*999)
                  ENDIF
                CASE(BASIS_SIMPLEX_TYPE)                  
                  CALL FLAG_ERROR("Regular meshes with simplex basis types is not implemented",ERR,ERROR,*999)
                CASE DEFAULT
                  CALL FLAG_ERROR("Basis type is either invalid or not implemented",ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Size of the mesh extent does not match the region number of dimensions",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Size of origin does not match the region number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Non rectangular cartesian coordinate systems are not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Region coordinate system is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_CREATE_REGULAR")
    RETURN
999 IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    CALL ERRORS("MESH_CREATE_REGULAR",ERR,ERROR)
    CALL EXITS("MESH_CREATE_REGULAR")
    RETURN 1   
  END SUBROUTINE MESH_CREATE_REGULAR
  
  !
  !================================================================================================================================
  !
  
END MODULE GENERATED_MESH_ROUTINES
