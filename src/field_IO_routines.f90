!> \file
!> $Id: parallel_IO.f90 1 2009-02-19 13:57:57Z cpb $
!> \author Heye Zhang
!> \brief ThiS module handles parallel Io. Using mpi2 and parall print function, 
!> formatted text and binary IO are supported in openCMISS  
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

!> Implements lists of Field IO operation
MODULE FIELD_IO_ROUTINES
  USE BASE_ROUTINES
  USE LISTS
  USE BASIS_ROUTINES
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE ISO_VARYING_STRING
  USE REGION_ROUTINES
  USE MACHINE_CONSTANTS
  USE KINDS
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE STRINGS
  USE TYPES
  USE CONSTANTS
  USE MPI
  USE CMISS_MPI
  USE INPUT_OUTPUT
  
  IMPLICIT NONE

  PRIVATE
  
  !Module parameters
  
  !> size of shape
  INTEGER(INTG), PARAMETER :: SHAPE_SIZE=3 

  !>Type for lable
  INTEGER(INTG), PARAMETER :: FIELD_IO_FIELD_LABEL=1  
  INTEGER(INTG), PARAMETER :: FIELD_IO_VARIABLE_LABEL=2  
  INTEGER(INTG), PARAMETER :: FIELD_IO_COMPONENT_LABEL=3  
  INTEGER(INTG), PARAMETER :: FIELD_IO_DERIVATIVE_LABEL=4

  !>Type of scale factor
  INTEGER(INTG), PARAMETER :: FIELD_IO_SCALE_FACTORS_NUMBER_TYPE=5
  INTEGER(INTG), PARAMETER :: FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE=6

  !Module types

  !>field variable compoment type pointer for IO
  TYPE MESH_ELEMENTS_TYPE_PTR_TYPE
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: PTR !< pointer field variable component
  END TYPE MESH_ELEMENTS_TYPE_PTR_TYPE

  !>field variable compoment type pointer for IO
  TYPE FIELD_VARIABLE_COMPONENT_PTR_TYPE
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: PTR !< pointer field variable component
  END TYPE FIELD_VARIABLE_COMPONENT_PTR_TYPE
  
  !>contains information for parallel IO, and it is nodal base
  TYPE FIELD_IO_INFO_SET
    !INTEGER(INTG) :: LEN_OF_NODAL_INFO !<how many bytes of this nodal information for IO
    LOGICAL :: SAME_HEADER !< determine whether we have same IO information as the previous one
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS !< number of components in the component array, COMPONENT(:)
    !attention: the pointers in COMPONENTS(:) point to those nodal components which are in the same local domain in current implementation
    !it may be replaced in the future implementation
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), ALLOCATABLE:: COMPONENTS(:) !<A array of pointers to those components of the node in this local domain       
  END TYPE FIELD_IO_INFO_SET  
  
  !>contains information for parallel IO, and it is nodal base
  TYPE FIELD_IO_NODAL_INFO_SET
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointe   r to the fields defined on the region.
    INTEGER(INTG) :: NUMBER_OF_NODES !<Number of nodes in this computional node for NODAL_INFO_SET
    !Interesting thing: pointer here, also means dymanically allocated attibute 
    INTEGER(INTG), ALLOCATABLE:: LIST_OF_GLOBAL_NUMBER(:) !<the list of global numbering in each domain  
    TYPE(FIELD_IO_INFO_SET), ALLOCATABLE:: NODAL_INFO_SET(:)  !<A list of nodal information for IO.      
  END TYPE FIELD_IO_NODAL_INFO_SET
  
  !>contains information for parallel IO, and it is nodal base
  TYPE FIELD_IO_ELEMENTALL_INFO_SET
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<A pointer to the fields defined on the region.
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS !<Number of nodes in this computional node for NODAL_INFO_SET
    !Interesting thing: pointer here, also means dymanically allocated attibute 
    INTEGER(INTG), ALLOCATABLE:: LIST_OF_GLOBAL_NUMBER(:) !<the list of global numbering in each domain  
    TYPE(FIELD_IO_INFO_SET), ALLOCATABLE:: ELEMENTAL_INFO_SET(:)  !<A list of nodal information for IO.      
  END TYPE FIELD_IO_ELEMENTALL_INFO_SET         
  
  !Module variables
  
  !Interfaces 
  PUBLIC :: FIELD_IO_NODES_EXPORT, FIELD_IO_ELEMENTS_EXPORT, FIELD_IO_FILEDS_IMPORT
 

CONTAINS  

  !
  !================================================================================================================================
  !  

  !>Get the derivative information               
  FUNCTION FIELD_IO_DERIVATIVE_INFO(LINE, ERR, ERROR)
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(IN) :: LINE !<Text info
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::FIELD_IO_DERIVATIVE_INFO
   
    CALL ENTERS("FIELD_IO_DERIVATIVE_INFO",ERR,ERROR,*999)    
   

    IF("d/ds1"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1
    ELSE IF("d2/ds1ds1"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S1       
    ELSE IF("d/ds2"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S2
    ELSE IF("d2/ds2ds2"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S2_S2
    ELSE IF("d/ds3"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S2
    ELSE IF("d2/ds3ds3"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S3
    ELSE IF("d2/ds3ds3"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S3_S3
    ELSE IF("d2/ds1ds3"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S3
    ELSE IF("d2/ds2ds3"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S2_S3 
    ELSE IF("d3/ds1ds2ds3"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S2_S3
    ELSE IF("d/ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S4
    ELSE IF("d2/ds4ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S4_S4
    ELSE IF("d2/ds1ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S4 
    ELSE IF("d2/ds2ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S2_S4
    ELSE IF("d2/ds3ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S3_S4
    ELSE IF("d3/ds1ds2ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S2_S4
    ELSE IF("d3/ds1ds3ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S3_S4
    ELSE IF("d3/ds2ds3ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S2_S3_S4
    ELSE IF("d4/ds1ds2ds3ds4"==LINE) THEN
       FIELD_IO_DERIVATIVE_INFO=PART_DERIV_S1_S2_S3_S4
    ELSE 
       FIELD_IO_DERIVATIVE_INFO=-1
       CALL FLAG_ERROR("Could not recognize derivatives from input string",ERR,ERROR,*999)
    ENDIF
             
    CALL EXITS("FIELD_IO_DERIVATIVE_INFO")
    RETURN
999 CALL ERRORS("FIELD_IO_DERIVATIVE_INFO",ERR,ERROR)
    CALL EXITS("FIELD_IO_DERIVATIVE_INFO")
  END FUNCTION FIELD_IO_DERIVATIVE_INFO    

  !
  !================================================================================================================================
  !  
  
  !>Create decompsition
  SUBROUTINE FIELD_IO_CREATE_FIELDS(NAME, REGION, DECOMPOSITION, FIELD_VALUES_SET_TYPE, NUMBER_OF_FIELDS, &
    &USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER, MESH_COMPONENTS_OF_FIELD_COMPONENTS, COMPONENTS_IN_FIELDS, NUMBER_OF_EXNODE_FILES, &
    &MASTER_COMPUTATIONAL_NUMBER, my_computational_node_number, FIELD_SCALING_TYPE, ERR, ERROR, *)
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(IN) :: NAME
    TYPE(REGION_TYPE), POINTER :: REGION !<region
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !< decompistion
    INTEGER(INTG), INTENT(IN) :: FIELD_VALUES_SET_TYPE !<
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_FIELDS !<!< number of fields
    INTEGER(INTG), INTENT(IN) :: USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(:)
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENTS_OF_FIELD_COMPONENTS(:)    
    INTEGER(INTG), INTENT(IN) :: COMPONENTS_IN_FIELDS(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_EXNODE_FILES
    INTEGER(INTG), INTENT(IN) :: MASTER_COMPUTATIONAL_NUMBER
    INTEGER(INTG), INTENT(IN) :: my_computational_node_number
    INTEGER(INTG), INTENT(IN) :: FIELD_SCALING_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<field
    TYPE(VARYING_STRING) :: FILE_NAME, FILE_STATUS, LINE, LINE1
    TYPE(VARYING_STRING) :: CMISS_KEYWORD_FIELDS, CMISS_KEYWORD_NODE, CMISS_KEYWORD_COMPONENTS
    TYPE(VARYING_STRING) :: CMISS_KEYWORD_VALUE_INDEX, CMISS_KEYWORD_DERIVATIVE
    INTEGER(INTG), ALLOCATABLE :: tmp_pointer(:), LIST_DEV(:), LIST_DEV_POS(:)
    INTEGER(INTG) :: FILE_ID
    INTEGER(INTG) :: NUMBER_FIELDS
    INTEGER(INTG) :: NODAL_USER_NUMBER
    INTEGER(INTG) :: MPI_IERROR
    INTEGER(INTG) :: idx_comp, idx_comp1, pos, idx_field, idx_exnodes
    INTEGER(INTG) :: idx_variable, idx_variable1, idx_dev, idx_dev1, total_number_of_comps, total_number_of_devs, number_of_devs
    INTEGER(INTG) :: number_of_coms
    REAL(DP), ALLOCATABLE :: LIST_DEV_VALUE(:)
    LOGICAL :: SECTION_START, FILE_END, NODE_SECTION
    
     
    CALL ENTERS("FIELD_IO_CREATE_FIELDS",ERR,ERROR,*999)   		
	
    IF(.NOT.ASSOCIATED(DECOMPOSITION)) THEN
      CALL FLAG_ERROR("decomposition is NOT associated before importing data",ERR,ERROR,*999)
      GOTO 999                     
    ENDIF

    IF(.NOT.ASSOCIATED(REGION)) THEN
      CALL FLAG_ERROR("region is NOT associated before importing data",ERR,ERROR,*999)
      GOTO 999                     
    ENDIF
    
    CMISS_KEYWORD_FIELDS="#Fields="
    CMISS_KEYWORD_COMPONENTS="#Components="
    CMISS_KEYWORD_VALUE_INDEX="Value index="
    CMISS_KEYWORD_DERIVATIVE="#Derivatives="
    CMISS_KEYWORD_NODE="Node:"
    
    idx_variable1=0
    idx_comp1=0
    DO idx_field=1,NUMBER_OF_FIELDS
       IF(ASSOCIATED(FIELD)) NULLIFY(FIELD)
       !Start to create a default (geometric) field on the region
       CALL FIELD_CREATE_START(idx_field,REGION,FIELD,ERR,ERROR,*999)       
       !always has one field variable in one field during reading
       CALL FIELD_NUMBER_OF_VARIABLES_SET(FIELD,1,ERR,ERROR,*999)
       !Set the decomposition to use
       CALL FIELD_MESH_DECOMPOSITION_SET(FIELD,DECOMPOSITION,ERR,ERROR,*999)       
       !Set the number of components for this field
       CALL FIELD_NUMBER_OF_COMPONENTS_SET(FIELD,COMPONENTS_IN_FIELDS(idx_field),ERR,ERROR,*999)       
       DO idx_comp=1, COMPONENTS_IN_FIELDS(idx_field)
          idx_comp1=idx_comp1+1
          !Set the domain to be used by the field components
          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(FIELD,1,idx_comp,MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp1), ERR,ERROR,*999)
       ENDDO   
       !Set the scaling factor
       CALL FIELD_SCALING_TYPE_SET(FIELD, FIELD_SCALING_TYPE, ERR, ERROR, *999)
       !Finish creating the field
       CALL FIELD_CREATE_FINISH(REGION,FIELD,ERR,ERROR,*999)
    ENDDO   
    
    DO WHILE(idx_exnodes<NUMBER_OF_EXNODE_FILES)              
       
       !goto the start of mesh part
       IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN
          
          IF(FILE_END==.TRUE.) THEN             
             FILE_ID=1030+idx_exnodes
             !checking the next file
             FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exnodes,"*",ERR,ERROR))//".exnode"
             CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
             FILE_END=.FALSE.
             SECTION_START=.FALSE.
             NODE_SECTION=.FALSE.
             idx_exnodes=idx_exnodes+1
          ENDIF
                    
          IF(FILE_END==.FALSE..AND.SECTION_START==.FALSE.)  THEN
             !find a new header
             DO WHILE(VERIFY(CMISS_KEYWORD_FIELDS,LINE)/=0)
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)                
             ENDDO
             SECTION_START=.TRUE.
          ENDIF
                    
          !have not touched the end
          IF(FILE_END==.FALSE..AND.SECTION_START==.TRUE..AND.NODE_SECTION==.FALSE.) THEN
             pos=INDEX(LINE,CMISS_KEYWORD_FIELDS)
             LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_FIELDS)-1)
             !number_of_fields=STRING_TO_INTEGER(LINE, ERR, ERROR)             
             total_number_of_comps=0
             total_number_of_devs=0
             idx_comp1=0
             idx_dev1=0
             DO idx_field=1, NUMBER_OF_FIELDS
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR, ERROR,*999)
                pos=INDEX(LINE,CMISS_KEYWORD_COMPONENTS)
                LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_COMPONENTS)-1)
                number_of_coms=STRING_TO_INTEGER(LINE, ERR, ERROR) 
                total_number_of_comps=total_number_of_comps+number_of_coms          
                
                IF(ALLOCATED(LIST_DEV_POS)) THEN                   
                   ALLOCATE(tmp_pointer(total_number_of_comps-number_of_coms),STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary memory for nodal derivative index in master node", &
                      &ERR,ERROR,*999)
                   tmp_pointer(:)=LIST_DEV_POS(:)
                   DEALLOCATE(LIST_DEV_POS)
                   ALLOCATE(LIST_DEV_POS(total_number_of_comps),STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary memory for nodal derivative index in master node", &
                      &ERR,ERROR,*999)
                   LIST_DEV_POS(1:total_number_of_comps-number_of_coms)=tmp_pointer(:)
                   DEALLOCATE(tmp_pointer)            
                ELSE
                   ALLOCATE(LIST_DEV_POS(total_number_of_comps),STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate memory for nodal derivative index",ERR,ERROR,*999)
                ENDIF
                
                DO idx_comp=1, number_of_coms                   
                   CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR, ERROR,*999)
                   pos=INDEX(LINE,".")
                   LINE=REMOVE(LINE,1, pos+1)
                   pos=INDEX(LINE,CMISS_KEYWORD_VALUE_INDEX)
                   LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_VALUE_INDEX)-1)
                   pos=INDEX(LINE,",")
                   LINE1=EXTRACT(LINE,1,pos-1)
                   idx_comp1=idx_comp1+1	
                   LIST_DEV_POS(idx_comp1)=STRING_TO_INTEGER(LINE1, ERR, ERROR)

                   pos=INDEX(LINE,CMISS_KEYWORD_DERIVATIVE)
                   LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_DERIVATIVE)-1)
                   pos=INDEX(LINE,"(")
                   LINE1=EXTRACT(LINE,1,pos-1)
                   number_of_devs=STRING_TO_INTEGER(LINE, ERR, ERROR)+1
                   total_number_of_devs=total_number_of_devs+number_of_devs  
                  
                   IF(ALLOCATED(LIST_DEV)) THEN                   
                      ALLOCATE(tmp_pointer(total_number_of_devs-number_of_devs),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary memory for nodal derivative index in master node", &
                         &ERR,ERROR,*999)
                      tmp_pointer(:)=LIST_DEV_POS(:)
                      DEALLOCATE(LIST_DEV_POS)
                      ALLOCATE(LIST_DEV_POS(total_number_of_devs),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary memory for nodal derivative index in master node", &
                         &ERR,ERROR,*999)
                      LIST_DEV_POS(1:total_number_of_devs-number_of_devs)=tmp_pointer(:)
                      DEALLOCATE(tmp_pointer)            
                   ELSE
                      ALLOCATE(LIST_DEV(total_number_of_devs),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate memory for nodal derivative index",ERR,ERROR,*999)
                   ENDIF
                                      
                   pos=INDEX(LINE,"(")
                   LINE=REMOVE(LINE,1, pos)
                   pos=INDEX(LINE,")")
                   LINE=REMOVE(LINE,pos, LEN(LINE))                                      
                   idx_dev1=idx_dev1+1
                   LIST_DEV(idx_dev1)=NO_PART_DERIV
                   DO idx_dev=2, number_of_devs-1
                      idx_dev1=idx_dev1+1
                      pos=INDEX(LINE,",")
                      LINE1=EXTRACT(LINE, 1, pos-1)
                      LINE=REMOVE(LINE, 1, pos)
                      LIST_DEV(idx_dev1)=FIELD_IO_DERIVATIVE_INFO(LINE1, ERR,ERROR)
                   ENDDO
                   idx_dev1=idx_dev1+1
                   LIST_DEV(idx_dev1)=FIELD_IO_DERIVATIVE_INFO(LINE, ERR,ERROR)                                        
                ENDDO !idx_comp
                NODE_SECTION=.TRUE.              
             ENDDO !idx_field
          ENDIF  !FILE_END==.FALSE..AND.SECTION_START=.TRUE..AND.NODE_SECTION=.FALSE.                  
       ENDIF !MASTER_COMPUTATIONAL_NUMBER

       !broadcasting total_number_of_comps
       CALL MPI_BCAST(total_number_of_comps,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                           
       !broadcasting total_number_of_devs
       CALL MPI_BCAST(total_number_of_devs,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)  
       
       IF(MASTER_COMPUTATIONAL_NUMBER/=my_computational_node_number) THEN
          IF(ALLOCATED(LIST_DEV_POS)) DEALLOCATE(LIST_DEV_POS)
          ALLOCATE(LIST_DEV_POS(total_number_of_comps),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate memory for nodal derivative index in non-master node",ERR, ERROR,*999)
          
          IF(ALLOCATED(LIST_DEV)) DEALLOCATE(LIST_DEV)
          ALLOCATE(LIST_DEV(total_number_of_devs),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate memory for nodal derivative index in non-master node",ERR, ERROR,*999)   
          
          IF(ALLOCATED(LIST_DEV_VALUE)) DEALLOCATE(LIST_DEV_VALUE)
          ALLOCATE(LIST_DEV_VALUE(total_number_of_devs),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate memory for nodal derivative index in non-master node",ERR, ERROR,*999)                    
       ENDIF

       !broadcasting total_number_of_comps
       CALL MPI_BCAST(LIST_DEV_POS,total_number_of_comps,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                           
       !broadcasting total_number_of_devs
       CALL MPI_BCAST(LIST_DEV,total_number_of_devs,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                

       !goto the start of mesh part
       IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN
                              
          !have not touched the end
          IF(FILE_END==.FALSE..AND.SECTION_START==.TRUE..AND.NODE_SECTION==.TRUE.) THEN
             CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR, ERROR,*999)
             IF(VERIFY(CMISS_KEYWORD_NODE, LINE)==0) THEN
                pos=INDEX(LINE,CMISS_KEYWORD_NODE)
                LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_NODE)-1)
                NODAL_USER_NUMBER=STRING_TO_INTEGER(LINE, ERR, ERROR)
                idx_comp1=1
                DO idx_comp=1, number_of_coms-1
                   CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR, ERROR,*999)
                   CALL STRING_TO_MUTI_REALS_VS(LINE, LIST_DEV_POS(idx_comp+1)-LIST_DEV_POS(idx_comp), LIST_DEV_VALUE, ERR,ERROR, *999)    
                ENDDO
                CALL STRING_TO_MUTI_REALS_VS(LINE, total_number_of_devs-LIST_DEV_POS(idx_comp)+1, LIST_DEV_VALUE, ERR,ERROR, *999)                
             ELSE
                NODE_SECTION=.TRUE.
             ENDIF !(VERIFY(CMISS_KEYWORD_NODE , LINE)==0)
          ENDIF !FILE_END==.FALSE..AND.SECTION_START=.TRUE..AND.NODE_SECTION=.TRUE.
       ENDIF  !(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number)  
       
       !broadcasting total_number_of_devs
       CALL MPI_BCAST(LIST_DEV_VALUE,total_number_of_devs,MPI_REAL,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
       
       idx_comp1=0 
       idx_dev1=0  
       DO idx_field=1,NUMBER_OF_FIELDS
          IF(ASSOCIATED(FIELD)) NULLIFY(FIELD)
          FIELD=>REGION%FIELDS%FIELDS(idx_field)%PTR
          DO idx_comp=1, COMPONENTS_IN_FIELDS(idx_field)
             idx_comp1=idx_comp1+1
             DO idx_dev=1, LIST_DEV_POS(idx_comp1+1)-LIST_DEV_POS(idx_comp1)
                idx_dev1=idx_dev1+1
                !Set the domain to be used by the field components
                CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_VALUES_SET_TYPE, LIST_DEV(idx_dev1), & 
                   &USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(NODAL_USER_NUMBER), idx_comp, idx_variable, LIST_DEV_VALUE(idx_dev1),&
                   &ERR, ERROR, *999)                        
             ENDDO !idx_dev       
          ENDDO !idx_comp  
       ENDDO !idx_field      
    ENDDO !idx_exnodes<NUMBER_OF_EXELEM_FILES  

    DO idx_field=1,NUMBER_FIELDS
       IF(ASSOCIATED(FIELD)) NULLIFY(FIELD)
       FIELD=>REGION%FIELDS%FIELDS(idx_field)%PTR
       CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
       CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
    ENDDO
    
    IF(ALLOCATED(tmp_pointer)) DEALLOCATE(tmp_pointer)
    IF(ALLOCATED(LIST_DEV_VALUE)) DEALLOCATE(LIST_DEV_VALUE)
    IF(ALLOCATED(LIST_DEV)) DEALLOCATE(LIST_DEV)
    IF(ALLOCATED(LIST_DEV_POS)) DEALLOCATE(LIST_DEV_POS)
 	
    CALL EXITS("FIELD_IO_CREATE_FIELDS")
    RETURN
999 CALL ERRORS("FIELD_IO_CREATE_FIELDS",ERR,ERROR)
    CALL EXITS("FIELD_IO_CREATE_FIELDS")
  END SUBROUTINE FIELD_IO_CREATE_FIELDS    


  !
  !================================================================================================================================
  !  
  
  !>Create decompition
  SUBROUTINE FIELD_IO_CREATE_DECOMPISTION(DECOMPOSITION, DECOMPOSITION_USER_NUMBER, DECOMPOSITION_METHOD, MESH, NUMBER_OF_DOMAINS, &
    &ERR, ERROR, *)
    !Argument variables   
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !< decomposition tye
    INTEGER(INTG), INTENT(IN) :: DECOMPOSITION_USER_NUMBER !<user number for decompistion
    INTEGER(INTG), INTENT(IN) :: DECOMPOSITION_METHOD !<decompistion type
    TYPE(MESH_TYPE), POINTER :: MESH !< mesh type
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !< number of domains 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    IF(.NOT.ASSOCIATED(MESH)) THEN
      CALL FLAG_ERROR("mesh is NOT associated before decomposing the mesh",ERR,ERROR,*999)
      GOTO 999                     
    ENDIF
    
    CALL ENTERS("FIELD_IO_CREATE_DECOMPISTION",ERR,ERROR,*999)   		
    !Create a decomposition
    CALL DECOMPOSITION_CREATE_START(DECOMPOSITION_USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_METHOD,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)
	                     
    CALL EXITS("FIELD_IO_CREATE_DECOMPISTION")
    RETURN
999 CALL ERRORS("FIELD_IO_CREATE_DECOMPISTION",ERR,ERROR)
    CALL EXITS("FIELD_IO_CREATE_DECOMPISTION")
  END SUBROUTINE FIELD_IO_CREATE_DECOMPISTION    


  !
  !================================================================================================================================
  !  
  
  !>Import fields from files into different computational nodes
  SUBROUTINE FIELD_IO_FILEDS_IMPORT(NAME, METHOD, REGION, MESH, MESH_USER_NUMBER, DECOMPOSITION, DECOMPOSITION_USER_NUMBER, &
    &DECOMPOSITION_METHOD, FIELD_VALUES_SET_TYPE, FIELD_SCALING_TYPE, ERR, ERROR, *)
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(IN) :: NAME !<name of input 
    TYPE(VARYING_STRING), INTENT(IN) :: METHOD !<method used for import
    TYPE(REGION_TYPE), POINTER :: REGION !<region
    TYPE(MESH_TYPE), POINTER :: MESH !<mesh type
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !<user number for mesh
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !< decompistion 
    INTEGER(INTG), INTENT(IN) :: DECOMPOSITION_USER_NUMBER !<user number for decompistion
    INTEGER(INTG), INTENT(IN) :: DECOMPOSITION_METHOD !<decompistion method
    INTEGER(INTG), INTENT(IN) :: FIELD_VALUES_SET_TYPE
    INTEGER(INTG), INTENT(IN) :: FIELD_SCALING_TYPE
    !TYPE(BASIS_FUNCTIONS_TYPE), POINTER :: BASES !< bases function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: my_computational_node_number !local process number    
    INTEGER(INTG) :: computational_node_numbers   !total process numbers
    INTEGER(INTG) :: MASTER_COMPUTATIONAL_NUMBER  !master computational number 
    INTEGER(INTG) :: NUMBER_OF_FIELDS
    INTEGER(INTG) :: NUMBER_OF_EXNODE_FILES
    INTEGER(INTG), ALLOCATABLE :: USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(:)
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENTS_OF_FIELD_COMPONENTS(:)
    INTEGER(INTG), ALLOCATABLE :: COMPONENTS_IN_FIELDS(:)
    CALL ENTERS("FIELD_IO_FILEDS_IMPORT",ERR,ERROR,*999)   		
       
    !Get the number of computational nodes
    computational_node_numbers=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !Get my computational node number
    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    
    MASTER_COMPUTATIONAL_NUMBER=0
    
    IF(METHOD=="FORTRAN") THEN
       CALL FIELD_IO_IMPORT_GLOBAL_MESH(NAME, REGION, MESH, MESH_USER_NUMBER, MASTER_COMPUTATIONAL_NUMBER, &
            & my_computational_node_number, USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER, MESH_COMPONENTS_OF_FIELD_COMPONENTS, &
            & COMPONENTS_IN_FIELDS, NUMBER_OF_FIELDS, NUMBER_OF_EXNODE_FILES, ERR, ERROR, *999)             
       
       CALL FIELD_IO_CREATE_DECOMPISTION(DECOMPOSITION, DECOMPOSITION_USER_NUMBER, DECOMPOSITION_METHOD, MESH, &
            &computational_node_numbers, ERR, ERROR, *999)
                 
       CALL FIELD_IO_CREATE_FIELDS(NAME, REGION, DECOMPOSITION, FIELD_VALUES_SET_TYPE, NUMBER_OF_FIELDS, &
            &USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER, MESH_COMPONENTS_OF_FIELD_COMPONENTS, COMPONENTS_IN_FIELDS, &
            &NUMBER_OF_EXNODE_FILES, MASTER_COMPUTATIONAL_NUMBER, my_computational_node_number, FIELD_SCALING_TYPE, ERR, ERROR, *999)
    ELSE IF(METHOD=="MPIIO") THEN
       CALL FLAG_ERROR("MPI IO has not been implemented",ERR,ERROR,*999)
    ELSE 
       CALL FLAG_ERROR("Unknown method!",ERR,ERROR,*999)   
    ENDIF     		
	
	IF(ALLOCATED(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER)) DEALLOCATE(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER)
	IF(ALLOCATED(MESH_COMPONENTS_OF_FIELD_COMPONENTS)) DEALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS)
	IF(ALLOCATED(COMPONENTS_IN_FIELDS)) DEALLOCATE(COMPONENTS_IN_FIELDS)
	                     
    CALL EXITS("FIELD_IO_FILEDS_IMPORT")
    RETURN
999 CALL ERRORS("FIELD_IO_FILEDS_IMPORT",ERR,ERROR)
    CALL EXITS("FIELD_IO_FILEDS_IMPORT")
  END SUBROUTINE FIELD_IO_FILEDS_IMPORT    

  !
  !================================================================================================================================
  !  
  
  !>Finding basis information 
  SUBROUTINE FIELD_IO_FILL_BASIS_INFO(INTERPOLATION_XI, LIST_STR, NUMBER_OF_COMPONENTS, ERR, ERROR, *)
    !Argument variables   
    INTEGER(INTG), INTENT(INOUT) :: INTERPOLATION_XI(:,:) !< xi interpolation type
    TYPE(VARYING_STRING), INTENT(INOUT) :: LIST_STR(:) !<label type
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS ! number of components
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LINE, LINE1
    INTEGER(INTG) :: idx_comp, pos
    INTEGER(INTG) :: num_interp, INTERPOLATION_TYPE
     
    CALL ENTERS("FIELD_IO_FILL_BASIS_INFO",ERR,ERROR,*999)   	
	
	DO idx_comp=1,NUMBER_OF_COMPONENTS
	   num_interp=0
	   LINE=LIST_STR(idx_comp)
       DO WHILE(VERIFY("*",LINE)==0)
          num_interp=num_interp+1 
          pos=INDEX(LINE,"*")
          LINE1=EXTRACT(LINE, 1, pos)
          LINE=REMOVE(LINE,1,pos)
          CALL FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE(INTERPOLATION_TYPE, LINE, ERR, ERROR, *999)
          INTERPOLATION_XI(idx_comp, num_interp)=INTERPOLATION_TYPE
       ENDDO
       num_interp=num_interp+1
       LINE1=EXTRACT(LINE, 1, pos)
       LINE=REMOVE(LINE,1,pos)
       CALL FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE(INTERPOLATION_TYPE, LINE, ERR, ERROR, *999)
       INTERPOLATION_XI(idx_comp, num_interp)=INTERPOLATION_TYPE	   
	ENDDO       
	                     
    CALL EXITS("FIELD_IO_FILL_BASIS_INFO")
    RETURN
999 CALL ERRORS("FIELD_IO_FILL_BASIS_INFO",ERR,ERROR)
    CALL EXITS("FIELD_IO_FILL_BASIS_INFO")
  END SUBROUTINE FIELD_IO_FILL_BASIS_INFO    


  !
  !================================================================================================================================
  !  
  
  !>Read the global mesh into one computational node first and then broadcasting to others nodes
  SUBROUTINE FIELD_IO_IMPORT_GLOBAL_MESH(NAME, REGION, MESH, MESH_USER_NUMBER, MASTER_COMPUTATIONAL_NUMBER, &
    & my_computational_node_number, USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER, MESH_COMPONENTS_OF_FIELD_COMPONENTS, &
    & COMPONENTS_IN_FIELDS, NUMBER_OF_FIELDS, NUMBER_OF_EXNODE_FILES, ERR, ERROR, *)
    !Argument variables       
    TYPE(VARYING_STRING), INTENT(IN):: NAME !< the name of elment file
    TYPE(MESH_TYPE), POINTER :: MESH !<mesh type
    TYPE(REGION_TYPE), POINTER :: REGION !<region
    INTEGER(INTG), INTENT(IN) :: MESH_USER_NUMBER !< user number of mesh      
    INTEGER(INTG), INTENT(IN) :: MASTER_COMPUTATIONAL_NUMBER 
    INTEGER(INTG), INTENT(IN) :: my_computational_node_number
    INTEGER(INTG), INTENT(INOUT), ALLOCATABLE :: USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(:)
    INTEGER(INTG), INTENT(INOUT), ALLOCATABLE :: MESH_COMPONENTS_OF_FIELD_COMPONENTS(:)
    INTEGER(INTG), INTENT(INOUT), ALLOCATABLE :: COMPONENTS_IN_FIELDS(:)
    INTEGER(INTG), INTENT(INOUT) :: NUMBER_OF_FIELDS
    INTEGER(INTG), INTENT(INOUT) :: NUMBER_OF_EXNODE_FILES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING), ALLOCATABLE :: LIST_STR(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(MESH_ELEMENTS_TYPE_PTR_TYPE), ALLOCATABLE :: ELEMENTS_PTR(:) 
    TYPE(VARYING_STRING) :: FILE_NAME, FILE_STATUS, LINE
    TYPE(VARYING_STRING) :: CMISS_KEYWORD_FIELD, CMISS_KEYWORD_ELEMENT, CMISS_KEYWORD_NODE, CMISS_KEYWORD_COMPONENT
    TYPE(VARYING_STRING) :: CMISS_KEYWORD_SHAPE, CMISS_KEYWORD_SCALE_FACTOR_SETS, CMISS_KEYWORD_NODES, CMISS_KEYWORD_SCALE_FACTORS
    INTEGER(INTG), PARAMETER :: NUMBER_NODAL_LINES=3, NUMBER_SCALING_FACTORS_IN_LINE=5
    INTEGER(INTG), ALLOCATABLE :: LIST_ELEMENT_NUMBER(:), LIST_ELEMENTAL_NODES(:), LIST_COMP_NODAL_INDEX(:,:)
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENT_LOOKUP(:,:), INTERPOLATION_XI(:,:), LIST_COMP_NODES(:)!LIST_FIELD_COMPONENTS(:),    
    INTEGER(INTG) :: FILE_ID, NUMBER_OF_EXELEM_FILES, NUMBER_OF_ELEMENTS, NUMBER_OF_NODES, NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: NUMBER_OF_MESH_COMPONENTS, NUMBER_OF_COMPONENTS, NUMBER_SCALING_FACTOR_LINES
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER
    INTEGER(INTG) :: MPI_IERROR
    INTEGER(INTG) :: SHAPE_INDEX(SHAPE_SIZE)
    INTEGER(INTG) :: idx_comp, idx_comp1, pos, idx_node, idx_node1, idx_field, idx_elem, idx_exnodes, idx_exelems, number_of_comp
    INTEGER(INTG) :: idx_basis, number_of_node, number_of_scalesets, idx_scl, idx_mesh_comp, current_mesh_comp, num_scl, num_scl_line
    LOGICAL :: FILE_EXIST, START_OF_ELEMENT_SECTION, FIELD_SECTION, SECTION_START, FILE_END
    
    CALL ENTERS("FIELD_IO_IMPORT_GLOBAL_MESH",ERR,ERROR,*999)    

    !checking region pointer
    IF(.NOT.ASSOCIATED(REGION)) THEN
       CALL FLAG_ERROR("region is not associated",ERR,ERROR,*999)
       GOTO 999         
    ENDIF

    !checking mesh pointer
    IF(ASSOCIATED(MESH)) THEN
       CALL FLAG_ERROR("mesh is associated, pls release the memory first",ERR,ERROR,*999)
       GOTO 999                 
    ENDIF
    
    IF(REGION%REGION_FINISHED==.FALSE.) THEN
       CALL FLAG_ERROR("region is not finished",ERR,ERROR,*999)
       GOTO 999         
    ENDIF
    
    !IF(ASSOCIATED(BASES%BASES)) THEN
    !   CALL FLAG_ERROR("bases are associated, pls release the memory first",ERR,ERROR,*999)
    !   GOTO 999
    !ENDIF
    !BASES%NUMBER_BASIS_FUNCTIONS=0            
    NUMBER_OF_DIMENSIONS=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
    FILE_STATUS="OLD"    
    CMISS_KEYWORD_SHAPE="Shape.  Dimension="//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))
    CMISS_KEYWORD_ELEMENT="Element:"
    CMISS_KEYWORD_COMPONENT="#Components="
    CMISS_KEYWORD_NODE="Node:"
    CMISS_KEYWORD_NODES="#Nodes="
    CMISS_KEYWORD_FIELD="#Fields="
    CMISS_KEYWORD_SCALE_FACTOR_SETS="#Scale factor sets="
    CMISS_KEYWORD_SCALE_FACTORS="#Scale factors="
    
    CALL MESH_CREATE_START(MESH_USER_NUMBER,REGION,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
    
    !calculate the number of elements, number of fields and number of field components
    IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN  
      
      !the file name has to start from zero in an ascended order without break      
      idx_exelems=0
      idx_elem=0
      NUMBER_OF_COMPONENTS=0
      FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exelems,"*",ERR,ERROR))//".exelem"
      INQUIRE(FILE=CHAR(FILE_NAME), EXIST=FILE_EXIST)
      IF(FILE_EXIST==.FALSE.) THEN
         CALL FLAG_ERROR("no file can be found, pls check again",ERR,ERROR,*999)
         GOTO 999                       
      ENDIF 
      DO WHILE(FILE_EXIST==.TRUE.)
         
         FILE_ID=1030+idx_exelems
         CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)         
         START_OF_ELEMENT_SECTION=.FALSE.
         FIELD_SECTION=.FALSE.
         
         CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR,*999) 
         DO WHILE(FILE_END==.FALSE.)                                                                                                    
            
            !check the beginning of element section
            IF(START_OF_ELEMENT_SECTION==.FALSE..AND.(VERIFY(CMISS_KEYWORD_SHAPE,LINE)==0)) THEN
               START_OF_ELEMENT_SECTION=.TRUE.               
            ENDIF   
            
            !count how many elements            
            IF(START_OF_ELEMENT_SECTION==.TRUE..AND.(VERIFY(CMISS_KEYWORD_ELEMENT,LINE)==0)) idx_elem=idx_elem+1
            
            !check whether they have same numbers of fields
            IF(START_OF_ELEMENT_SECTION==.TRUE..AND.(VERIFY(CMISS_KEYWORD_FIELD,LINE)==0)) THEN
               idx_field=0
               idx_comp=0
               FIELD_SECTION=.TRUE.
               pos=INDEX(LINE,CMISS_KEYWORD_FIELD)               
               LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_FIELD)-1)               
               IF(idx_exelems==0) THEN
                  NUMBER_OF_FIELDS=STRING_TO_INTEGER(LINE, ERR,ERROR)
               ELSE
                  IF(NUMBER_OF_FIELDS/=STRING_TO_INTEGER(LINE, ERR,ERROR)) THEN
                     CALL FLAG_ERROR("find different number of fields in the exelem files",ERR,ERROR,*999)
                     GOTO 999
                  ENDIF
               ENDIF

               IF(.NOT.ALLOCATED(COMPONENTS_IN_FIELDS)) THEN
                  ALLOCATE(COMPONENTS_IN_FIELDS(NUMBER_OF_FIELDS), STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("can not allocate the momery for outputing components in field",ERR,ERROR,*999)
                  COMPONENTS_IN_FIELDS(:)=0
               ENDIF   
            ENDIF !START_OF_ELEMENT_SECTION==.TRUE..AND.VERIFY(CMISS_KEYWORD_FIELD,LINE)==0  
            
            !check whether they have same numbers of field components
            IF(FIELD_SECTION==.TRUE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.(VERIFY(CMISS_KEYWORD_COMPONENT,LINE)==0)) THEN               
               idx_field=idx_field+1
               pos=INDEX(LINE,CMISS_KEYWORD_COMPONENT)               
               LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_COMPONENT)-1)
               idx_comp1=STRING_TO_INTEGER(LINE, ERR,ERROR)
               idx_comp=idx_comp+idx_comp1
               IF(idx_field>=NUMBER_OF_FIELDS) THEN   
                  IF(idx_exelems==0) THEN
                     NUMBER_OF_COMPONENTS=idx_comp
                     COMPONENTS_IN_FIELDS(idx_field)=idx_comp
                  ELSE
                     IF(NUMBER_OF_COMPONENTS/=idx_comp) THEN
                        CALL FLAG_ERROR("find different total number of components in the exelem files",ERR,ERROR,*999)
                        GOTO 999
                     ENDIF                     
                  ENDIF   
                  FIELD_SECTION=.FALSE.
               ENDIF
               IF(idx_exelems==0) THEN
                 COMPONENTS_IN_FIELDS(idx_field)=idx_comp1
               ELSE
                 IF(COMPONENTS_IN_FIELDS(idx_field)/=idx_comp1) THEN
                    CALL FLAG_ERROR("find different number of components in one field in the exelem files",ERR,ERROR,*999)
                    GOTO 999               
                 ENDIF
               ENDIF   
            ENDIF !FIELD_SECTION==.TRUE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.VERIFY(CMISS_KEYWORD_COMPONENT,LINE
            CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR,*999)                                       
         ENDDO !(FILE_END==.FALSE.)                  
                  
         CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
         !checking the next file
         idx_exelems=idx_exelems+1
         FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exelems,"*",ERR,ERROR))//".exelem"
         INQUIRE(FILE=CHAR(FILE_NAME), EXIST=FILE_EXIST)
      ENDDO!FILE_EXIST==.TRUE.
      !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Total number of exelment files = ",idx_exelems, ERR,ERROR,*999)
      NUMBER_OF_ELEMENTS=idx_elem
      NUMBER_OF_EXELEM_FILES=idx_exelems
    ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number           

        
    !broadcasting the number of components in each field
    CALL MPI_BCAST(NUMBER_OF_FIELDS,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
    IF(MASTER_COMPUTATIONAL_NUMBER/=my_computational_node_number) THEN 
       IF(ALLOCATED(COMPONENTS_IN_FIELDS)) DEALLOCATE(COMPONENTS_IN_FIELDS)
       ALLOCATE(COMPONENTS_IN_FIELDS(NUMBER_OF_FIELDS), STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("can not allocate the momery for outputing components in field",ERR,ERROR,*999)
       COMPONENTS_IN_FIELDS(:)=0
    ENDIF
    CALL MPI_BCAST(COMPONENTS_IN_FIELDS,NUMBER_OF_FIELDS,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)    
    !broadcasting the number of elements
    CALL MPI_BCAST(NUMBER_OF_ELEMENTS,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*999)

    !calculate the number of nodes
    IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN     
      !the file name has to start from zero in a ascended order without break    
      idx_exnodes=0
      idx_node=0
      FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exnodes,"*",ERR,ERROR))//".exnode"
      INQUIRE(FILE=CHAR(FILE_NAME), EXIST=FILE_EXIST)
      IF(FILE_EXIST==.FALSE.) THEN
         CALL FLAG_ERROR("no file can be found, pls check again",ERR,ERROR,*999)
         GOTO 999                       
      ENDIF 
      DO WHILE(FILE_EXIST==.TRUE.)      
         FILE_ID=1030+idx_exnodes
         CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
         FILE_END=.FALSE.
         
         DO WHILE(FILE_END==.FALSE.)    
            CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
            IF(FILE_END==.FALSE..AND.VERIFY(CMISS_KEYWORD_NODE,LINE)==0) idx_node=idx_node+1            
         ENDDO !(FILE_END==.FALSE.)                 
         
         CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)         
         !checking the next file
         idx_exnodes=idx_exnodes+1
         FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exnodes,"*",ERR,ERROR))//".exnode"
         INQUIRE(FILE=CHAR(FILE_NAME), EXIST=FILE_EXIST)
      ENDDO!FILE_EXIST==.TRUE.
      !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Total number of exnode files = ",idx_exnodes, ERR,ERROR,*999)
      NUMBER_OF_NODES=idx_node
      NUMBER_OF_EXNODE_FILES=idx_exnodes
    ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number           
        
    CALL MPI_BCAST(NUMBER_OF_EXNODE_FILES,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
    !broadcasting the number of nodes
    CALL MPI_BCAST(NUMBER_OF_NODES,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
    CALL NODES_CREATE_START(NUMBER_OF_NODES,REGION,NODES,ERR,ERROR,*999)
    
    !collect the nodal numberings (nodal labels) to change the nodal user number by reading exnode files
    IF(ALLOCATED(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER)) DEALLOCATE(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER)    
    ALLOCATE(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(NUMBER_OF_NODES),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of nodal number",ERR,ERROR,*999)
    IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN    
      !the file name has to start from zero in a ascended order without break
      idx_node=1
      DO idx_exnodes=0, NUMBER_OF_EXNODE_FILES-1
         FILE_ID=1030+idx_exnodes
         FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exnodes,"*",ERR,ERROR))//".exnode"
         CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)                  
         FILE_END=.FALSE.
         DO WHILE(FILE_END==.FALSE.)    
            CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
            IF(FILE_END==.FALSE..AND.VERIFY(CMISS_KEYWORD_NODE,LINE)==0) THEN
               pos=INDEX(LINE,CMISS_KEYWORD_NODE)               
               LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_NODE)-1)
               USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(idx_node)=STRING_TO_INTEGER(LINE, ERR, ERROR)
               idx_node=idx_node+1
            ENDIF!VERIFY(CMISS_KEYWORD,LINE)==0
         ENDDO !(FILE_END==.FALSE.)                          
         CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
      ENDDO!FILE_EXIST==.TRUE.
      CALL LIST_SORT(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER, ERR, ERROR, *999)
    ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number           

    !broadcast the nodal numberings (nodal labels)
    CALL MPI_BCAST(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER,NUMBER_OF_NODES,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)         
    DO idx_node=1, NUMBER_OF_NODES
       IF(idx_node/=USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(idx_node)) CALL NODE_NUMBER_SET(idx_node, USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(idx_node), NODES, ERR,ERROR,*999)
    ENDDO        
    CALL NODES_CREATE_FINISH(REGION,ERR,ERROR,*999)
        
    !IF(ALLOCATED(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER)) DEALLOCATE(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER)
    !ALLOCATE(USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(NUMBER_OF_NODES), STAT=ERR)
    !IF(ERR/=0) CALL FLAG_ERROR("can not allocate nodal number mapping for output",ERR,ERROR,*999)
    !USER_NODAL_NUMBER_MAP_GLOBAL_NODAL_NUMBER(:)=LIST_NODAL_NUMBER(:)
    !IF(ALLOCATED(LIST_NODAL_NUMBER)) DEALLOCATE(LIST_NODAL_NUMBER)

    !calculate the number of mesh components
    IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN      
      
      !MESH_COMPONENT_LOOKUP is used to store the difference between field components in term of basis property.
      IF(ALLOCATED(MESH_COMPONENT_LOOKUP)) DEALLOCATE(MESH_COMPONENT_LOOKUP)
      ALLOCATE(MESH_COMPONENT_LOOKUP(NUMBER_OF_COMPONENTS,NUMBER_OF_COMPONENTS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of mesh components",ERR,ERROR,*999)
      
      IF(ALLOCATED(LIST_STR)) DEALLOCATE(LIST_STR)
      ALLOCATE(LIST_STR(NUMBER_OF_COMPONENTS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of str",ERR,ERROR,*999)
      !initialize MESH_COMPONENT_LOOKUP and assume each field component has the same mesh component
      MESH_COMPONENT_LOOKUP(:,:)=0
      DO idx_comp=1,NUMBER_OF_COMPONENTS
         MESH_COMPONENT_LOOKUP(idx_comp,idx_comp)=1
      ENDDO
      idx_exelems=0
      
      !checking field component's mesh component by checking the basis           
      DO WHILE(idx_exelems<NUMBER_OF_EXELEM_FILES)
         
         FILE_ID=1030+idx_exelems
         !checking the next file
         FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exelems,"*",ERR,ERROR))//".exelem"
         CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)         
         FIELD_SECTION=.FALSE.    
         FILE_END=.FALSE. 
          
         DO WHILE(FILE_END==.FALSE.)    
            CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
            
            !check the beginning of element section
            IF(START_OF_ELEMENT_SECTION==.FALSE..AND.(VERIFY(CMISS_KEYWORD_SHAPE,LINE)==0)) THEN
               START_OF_ELEMENT_SECTION=.TRUE.               
            ENDIF   
                       
            !check whether it is a new header for another group of elements
            IF(START_OF_ELEMENT_SECTION==.TRUE..AND.(VERIFY(CMISS_KEYWORD_FIELD,LINE)==0)) THEN
               
               !collect header information
               pos=INDEX(LINE,CMISS_KEYWORD_FIELD)               
               LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_FIELD)-1)
               NUMBER_OF_FIELDS=STRING_TO_INTEGER(LINE, ERR,ERROR)
               idx_comp=0
               DO idx_field=1,NUMBER_OF_FIELDS
                  CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                  pos=INDEX(LINE,CMISS_KEYWORD_COMPONENT)               
                  LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_COMPONENT)-1)
                  number_of_comp=STRING_TO_INTEGER(LINE, ERR,ERROR)               
                  DO idx_comp1=1,number_of_comp
                     idx_comp=idx_comp+1
                     CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                     pos=INDEX(LINE,".")
                     LINE=REMOVE(LINE,1, pos)
                     LINE=TRIM(ADJUSTL(LINE))
                     pos=INDEX(LINE,",")                     
                     LINE=REMOVE(LINE,pos, LEN(LINE))
                     LIST_STR(idx_comp)= LINE
                     CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                     pos=INDEX(LINE, CMISS_KEYWORD_NODES)
                     LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_NODES)-1)
                     number_of_node=STRING_TO_INTEGER(LINE, ERR,ERROR)
                     DO idx_node1=1, number_of_node*NUMBER_NODAL_LINES
                        CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                     ENDDO !idx_node1                
                  ENDDO !idx_comp1
               ENDDO !idx_field
               
               !compare the bases. Since the geometrical topology is same, if the bases are the same in the same topology, 
               !they should be in the same mesh component
               IF(SUM(MESH_COMPONENT_LOOKUP)<NUMBER_OF_COMPONENTS*NUMBER_OF_COMPONENTS) THEN
                  DO idx_comp=1, NUMBER_OF_COMPONENTS
                     DO idx_comp1=idx_comp+1, NUMBER_OF_COMPONENTS 
                        IF(MESH_COMPONENT_LOOKUP(idx_comp1,idx_comp)==0) THEN
                           IF(LIST_STR(idx_comp1)/=LIST_STR(idx_comp)) THEN
                              MESH_COMPONENT_LOOKUP(idx_comp1,idx_comp)=1
                              MESH_COMPONENT_LOOKUP(idx_comp,idx_comp1)=1
                           ENDIF
                        ENDIF
                     ENDDO !idx_comp1   
                  ENDDO !idx_comp
               ELSE
                  idx_exelems=NUMBER_OF_EXELEM_FILES!jump out of the loop
               ENDIF    
            ENDIF !START_OF_ELEMENT_SECTION==.TRUE..AND.VERIFY(CMISS_KEYWORD_FIELD,LINE)==0  
            
            !!check whether they have same numbers of field components
            !IF(FIELD_SECTION==.TRUE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.(VERIFY(CMISS_KEYWORD_COMPONENT,LINE)==0)) THEN                             
            !   IF(idx_field>=NUMBER_OF_FIELDS) THEN   
            !      FIELD_SECTION=.FALSE.
            !   ENDIF   
            !ENDIF !FIELD_SECTION==.TRUE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.VERIFY(CMISS_KEYWORD_COMPONENT,LINE             
         ENDDO !(FILE_END==.FALSE.)                                    
         CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
         idx_exelems=idx_exelems+1
      ENDDO!idx_exelems=0, NUMBER_OF_EXELEM_FILES-1                            
      
      !calculate the number of mesh components
      IF(ALLOCATED(MESH_COMPONENTS_OF_FIELD_COMPONENTS)) DEALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS)
      ALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS(NUMBER_OF_COMPONENTS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of field components",ERR,ERROR,*999)
      DO idx_comp=1, NUMBER_OF_COMPONENTS
         MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp)=idx_comp
      ENDDO
      DO idx_comp=1, NUMBER_OF_COMPONENTS
         DO idx_comp1=idx_comp+1, NUMBER_OF_COMPONENTS
            IF(MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp)==idx_comp) THEN
               IF(MESH_COMPONENT_LOOKUP(idx_comp1,idx_comp)==0) MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp1)=idx_comp
            ENDIF   
         ENDDO !idx_comp1
      ENDDO !idx_comp
      IF(ALLOCATED(MESH_COMPONENT_LOOKUP)) DEALLOCATE(MESH_COMPONENT_LOOKUP)      
      NUMBER_OF_MESH_COMPONENTS=0
      idx_comp1=0          
      DO idx_comp=1,NUMBER_OF_COMPONENTS
         IF(MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp)==idx_comp) THEN
            idx_comp1=idx_comp1+1
            MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp)=idx_comp1
            NUMBER_OF_MESH_COMPONENTS=NUMBER_OF_MESH_COMPONENTS+1
         ENDIF
      ENDDO
      !IF(ALLOCATED(MESH_COMPONENTS_OF_FIELD_COMPONENTS)) DEALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS)
      !ALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS(NUMBER_OF_COMPONENTS),STAT=ERR)
      !IF(ERR/=0) CALL FLAG_ERROR("can not allocate mesh components of field components for output",ERR,ERROR,*999)
      !MESH_COMPONENTS_OF_FIELD_COMPONENTS(:)=LIST_FIELD_COMPONENTS(:)
      !ALLOCATE(LIST_MESH_COMPONENTS(NUMBER_OF_MESH_COMPONENTS),STAT=ERR)
      !IF(ERR/=0) CALL FLAG_ERROR("ALLOCATEDcan not allocate list of mesh components",ERR,ERROR,*999)
      !idx_comp1=0
      !DO idx_comp=1,NUMBER_OF_COMPONENTS        
      !   IF(LIST_FIELD_COMPONENTS(idx_comp)==idx_comp) THEN
      !   	   idx_comp1=idx_comp1+1
      !      LIST_MESH_COMPONENTS(idx_comp1)=idx_comp
      !   ENDIF   
      !ENDDO      
    ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number           

    !broadcasting the number of mesh components        
    CALL MPI_BCAST(NUMBER_OF_COMPONENTS,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
    CALL MPI_BCAST(NUMBER_OF_MESH_COMPONENTS,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
    CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_MESH_COMPONENTS,ERR,ERROR,*999)
    IF(ALLOCATED(ELEMENTS_PTR)) DEALLOCATE(ELEMENTS_PTR)     
    ALLOCATE(ELEMENTS_PTR(NUMBER_OF_MESH_COMPONENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of mesh element pointers",ERR,ERROR,*999)
    
    IF(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS<=0)  THEN
       CALL BASIS_CREATE_START(1,BASIS,ERR,ERROR,*999)
       CALL BASIS_NUMBER_OF_XI_SET(BASIS,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
       CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
    ENDIF   
    
    DO idx_comp=1, NUMBER_OF_MESH_COMPONENTS
       CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,idx_comp,BASIS_FUNCTIONS%BASES(1)%PTR,ELEMENTS_PTR(idx_comp)%PTR,ERR,ERROR,*999)
    ENDDO
    
    !Collect the elemental numberings (elemental labels)
    IF(ALLOCATED(LIST_ELEMENT_NUMBER)) DEALLOCATE(LIST_ELEMENT_NUMBER) 
    ALLOCATE(LIST_ELEMENT_NUMBER(NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of elemental number",ERR,ERROR,*999)		    
    IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN 
      
      !the file name has to start from zero in a ascended order without break
      idx_elem=1
      DO idx_exelems=0, NUMBER_OF_EXELEM_FILES-1
         
         FILE_ID=1030+idx_exelems
         FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exelems,"*",ERR,ERROR))//".exelem"      
         CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)         
         START_OF_ELEMENT_SECTION=.FALSE. 
         FILE_END=.FALSE.  

         
         DO WHILE(FILE_END==.FALSE.)    
            CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
            
            IF(START_OF_ELEMENT_SECTION==.FALSE..AND.(VERIFY(CMISS_KEYWORD_SHAPE,LINE)==0)) THEN
               START_OF_ELEMENT_SECTION=.TRUE.               
            ENDIF   

            IF(START_OF_ELEMENT_SECTION==.TRUE..AND.(VERIFY(CMISS_KEYWORD_ELEMENT,LINE)==0)) THEN
               pos=INDEX(LINE,CMISS_KEYWORD_ELEMENT)
               LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_ELEMENT)-1)
               SHAPE_INDEX(:)=0
               CALL STRING_TO_MUTI_INTEGERS_VS(LINE, SHAPE_SIZE, SHAPE_INDEX(:), ERR,ERROR, *999)
               IF(NUMBER_OF_DIMENSIONS==3) THEN
                  LIST_ELEMENT_NUMBER(idx_elem)=SHAPE_INDEX(1)
               ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
                  LIST_ELEMENT_NUMBER(idx_elem)=SHAPE_INDEX(2)
               ELSE IF(NUMBER_OF_DIMENSIONS==1) THEN
                  LIST_ELEMENT_NUMBER(idx_elem)=SHAPE_INDEX(3)
               ELSE
                  CALL FLAG_ERROR("Non recognized dimension size during reading elemental numbering",ERR,ERROR,*999)        
               ENDIF     
               idx_elem=idx_elem+1
            ENDIF                            
         ENDDO !(FILE_END==.FALSE.)                  
                  
         CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
      ENDDO!idx_exelems=0
      CALL LIST_SORT(LIST_ELEMENT_NUMBER, ERR, ERROR, *999)      
    ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number           
    
    !broadcast the list of elements for mapping gloabl numbers and user numbers (elemental labels)
    CALL MPI_BCAST(LIST_ELEMENT_NUMBER,NUMBER_OF_ELEMENTS,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                                      
    !change the mapping between global elemental numbering and user elemental numbering
        
    DO idx_elem=1,NUMBER_OF_ELEMENTS
       DO idx_comp=1, NUMBER_OF_MESH_COMPONENTS
          IF(idx_elem/=LIST_ELEMENT_NUMBER(idx_elem)) &            
            &CALL MESH_TOPOLOGY_ELEMENTS_NUMBER_SET(idx_elem,LIST_ELEMENT_NUMBER(idx_elem),ELEMENTS_PTR(idx_comp)%PTR,ERR,ERROR,*999)             
       ENDDO
    ENDDO      
    
    !creating topological information for each mesh component
    CALL MPI_BCAST(NUMBER_OF_EXELEM_FILES,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
    !ALLOCATE(LIST_BASES(NUMBER_OF_COMPONENTS),STAT=ERR)
    !IF(ERR/=0) CALL FLAG_ERROR("can not allocate list of bases",ERR,ERROR,*999)		    

    IF(ALLOCATED(LIST_COMP_NODES)) DEALLOCATE(LIST_COMP_NODES) 
    ALLOCATE(LIST_COMP_NODES(NUMBER_OF_COMPONENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of component nodal index ",ERR,ERROR,*999)             
    LIST_COMP_NODES(:)=0

    FILE_END=.TRUE.
    idx_exelems=-1        
    
    DO WHILE(idx_exelems<NUMBER_OF_EXELEM_FILES)              
       
       !IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN
       !   CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"FUCK U in file:",idx_exelems, ERR,ERROR,*999)
          !PRINT *, idx_exelems
       !ENDIF     
              
       CALL MPI_BCAST(FILE_END,1,MPI_LOGICAL,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                                              
       
       IF(FILE_END==.TRUE.) THEN
          idx_exelems=idx_exelems+1
          IF(idx_exelems>=NUMBER_OF_EXELEM_FILES) EXIT
       ENDIF
              
       !goto the start of mesh part
       IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN
                    
          !IF(FILE_END==.FALSE..AND.START_OF_ELEMENT_SECTION==.FALSE.) THEN
          !   !check the beginning of element section          
          !   DO WHILE(VERIFY(CMISS_KEYWORD_SHAPE,LINE)/=0)
          !      CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
          !   ENDDO
          !   START_OF_ELEMENT_SECTION=.TRUE.
          !ENDIF
          
          IF(FILE_END==.TRUE.) THEN
             FILE_ID=1030+idx_exelems
             !checking the next file             
             FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(idx_exelems,"*",ERR,ERROR))//".exelem"
             CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
             FILE_END=.FALSE.
             SECTION_START=.FALSE.
             START_OF_ELEMENT_SECTION=.FALSE.
          ENDIF        
          
          
          IF(FILE_END==.FALSE..AND.START_OF_ELEMENT_SECTION==.FALSE.) THEN !..AND.SECTION_START==.FALSE.) THEN    
             !find a new header
             DO WHILE(VERIFY(CMISS_KEYWORD_SCALE_FACTOR_SETS,LINE)/=0)
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)    
             ENDDO
             START_OF_ELEMENT_SECTION=.TRUE.
          ENDIF
                                        
          !have not touched the end
          IF(FILE_END==.FALSE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.SECTION_START==.FALSE.) THEN
             SECTION_START=.TRUE.
             !CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,LINE,ERR,ERROR,*999)
             !CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILE_NAME,ERR,ERROR,*999)
             
             pos=INDEX(LINE,CMISS_KEYWORD_SCALE_FACTOR_SETS)
             LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_SCALE_FACTOR_SETS)-1)
             number_of_scalesets=STRING_TO_INTEGER(LINE, ERR,ERROR)
             idx_mesh_comp=1
             !skip factors
             NUMBER_SCALING_FACTOR_LINES=0
             DO idx_scl=1,number_of_scalesets
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                pos=INDEX(LINE,CMISS_KEYWORD_SCALE_FACTORS)
                LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_SCALE_FACTORS)-1)
                num_scl=STRING_TO_INTEGER(LINE, ERR,ERROR)
                num_scl_line=num_scl/NUMBER_SCALING_FACTORS_IN_LINE
                IF(num_scl_line*NUMBER_SCALING_FACTORS_IN_LINE/=num_scl) num_scl_line=num_scl_line+1
                NUMBER_SCALING_FACTOR_LINES=NUMBER_SCALING_FACTOR_LINES+num_scl_line
             ENDDO
             !skip nodes line
             CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
             pos=INDEX(LINE,CMISS_KEYWORD_NODES)               
             LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_NODES)-1)
             number_of_node=STRING_TO_INTEGER(LINE, ERR,ERROR)
             
             IF(ALLOCATED(LIST_ELEMENTAL_NODES)) DEALLOCATE(LIST_ELEMENTAL_NODES)
             ALLOCATE(LIST_ELEMENTAL_NODES(number_of_node),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of elemental nodes",ERR,ERROR,*999)
             !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"number_of_node:",number_of_node,ERR,ERROR,*999)
             LIST_ELEMENTAL_NODES(:)=0  
             IF(ALLOCATED(LIST_COMP_NODAL_INDEX)) DEALLOCATE(LIST_COMP_NODAL_INDEX)       
             ALLOCATE(LIST_COMP_NODAL_INDEX(NUMBER_OF_COMPONENTS,number_of_node),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of component nodal index ",ERR,ERROR,*999)
             LIST_COMP_NODAL_INDEX(:,:)=0                          
             
             !read the header
             CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
             pos=INDEX(LINE,CMISS_KEYWORD_FIELD)               
             LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_FIELD)-1)
             NUMBER_OF_FIELDS=STRING_TO_INTEGER(LINE, ERR,ERROR)
             idx_comp=0
             DO idx_field=1,NUMBER_OF_FIELDS
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                pos=INDEX(LINE,CMISS_KEYWORD_COMPONENT)               
                LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_COMPONENT)-1)
                number_of_comp=STRING_TO_INTEGER(LINE, ERR,ERROR)               
                DO idx_comp1=1,number_of_comp
                   idx_comp=idx_comp+1
                   CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                   pos=INDEX(LINE,".")
                   LINE=REMOVE(LINE,1, pos)
                   LINE=TRIM(ADJUSTL(LINE))
                   pos=INDEX(LINE,",")                     
                   LINE=REMOVE(LINE,pos, LEN(LINE))
                   LIST_STR(idx_comp)= LINE                  
                   CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                   pos=INDEX(LINE, CMISS_KEYWORD_NODES)
                   LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_NODES)-1)
                   number_of_node=STRING_TO_INTEGER(LINE, ERR,ERROR)
                   LIST_COMP_NODES(idx_comp)=number_of_node
                   DO idx_node1=1, number_of_node
                      CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                      pos=INDEX(LINE, ".")
                      LINE=REMOVE(LINE,pos, LEN(LINE))
                      LIST_COMP_NODAL_INDEX(idx_comp,idx_node1)=STRING_TO_INTEGER(LINE, ERR,ERROR)
                      CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                      CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
                   ENDDO !idx_node1                
                ENDDO !idx_comp1                
             ENDDO !idx_field   
             
             IF(ALLOCATED(INTERPOLATION_XI)) DEALLOCATE(INTERPOLATION_XI)    
             ALLOCATE(INTERPOLATION_XI(NUMBER_OF_COMPONENTS,NUMBER_OF_DIMENSIONS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of interpolation types",ERR,ERROR,*999)

             CALL FIELD_IO_FILL_BASIS_INFO(INTERPOLATION_XI, LIST_STR, NUMBER_OF_COMPONENTS, ERR,ERROR,*999)             
             CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)             
             !CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,LINE,ERR,ERROR,*999)
          ENDIF  !FILE_END==.FALSE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.SECTION_START==.TRUE.                              
          
          IF(FILE_END==.FALSE..AND.START_OF_ELEMENT_SECTION==.TRUE..AND.SECTION_START==.TRUE.) THEN                          
             
             pos=INDEX(LINE,CMISS_KEYWORD_ELEMENT)
             LINE=REMOVE(LINE,1, pos+LEN_TRIM(CMISS_KEYWORD_ELEMENT)-1)
             CALL STRING_TO_MUTI_INTEGERS_VS(LINE, SHAPE_SIZE, SHAPE_INDEX(:), ERR,ERROR, *999)                          
             
             
             DO WHILE(VERIFY(CMISS_KEYWORD_NODE,LINE)/=0)
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)                
             ENDDO
             
             CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
             CALL STRING_TO_MUTI_INTEGERS_VS(LINE, number_of_node, LIST_ELEMENTAL_NODES, ERR, ERROR, *999)
             
             !skip scaling factors
             CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
             DO idx_scl=1, NUMBER_SCALING_FACTOR_LINES
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)                
             ENDDO
             
             IF(FILE_END==.FALSE.) THEN
                CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)             
                IF(VERIFY(CMISS_KEYWORD_SCALE_FACTOR_SETS,LINE)==0) SECTION_START=.TRUE.
             ENDIF   
             !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"FILE_END:",FILE_END,ERR,ERROR,*999)
          ENDIF                        
       ENDIF !MASTER_COMPUTATIONAL_NUMBER
       
       CALL MPI_BCAST(number_of_node,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                                      
       !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"SIZE LIST_ELEMENTAL_NODES:",SIZE(LIST_ELEMENTAL_NODES),ERR,ERROR,*999)
       
       IF(MASTER_COMPUTATIONAL_NUMBER/=my_computational_node_number) THEN
          IF(ALLOCATED(LIST_ELEMENTAL_NODES)) DEALLOCATE(LIST_ELEMENTAL_NODES)
          ALLOCATE(LIST_ELEMENTAL_NODES(number_of_node),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of elemental nodes",ERR,ERROR,*999)
          LIST_ELEMENTAL_NODES(:)=0  
          IF(ALLOCATED(LIST_COMP_NODAL_INDEX)) DEALLOCATE(LIST_COMP_NODAL_INDEX)       
          ALLOCATE(LIST_COMP_NODAL_INDEX(NUMBER_OF_COMPONENTS,number_of_node),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of component nodal index ",ERR,ERROR,*999)
          LIST_COMP_NODAL_INDEX(:,:)=0
          IF(ALLOCATED(INTERPOLATION_XI)) DEALLOCATE(INTERPOLATION_XI)       
          ALLOCATE(INTERPOLATION_XI(NUMBER_OF_COMPONENTS,NUMBER_OF_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of interpolation types",ERR,ERROR,*999)
          INTERPOLATION_XI(:,:)=0
          IF(ALLOCATED(MESH_COMPONENTS_OF_FIELD_COMPONENTS)) DEALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS)
          ALLOCATE(MESH_COMPONENTS_OF_FIELD_COMPONENTS(NUMBER_OF_COMPONENTS), STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of mesh components of field",ERR,ERROR,*999)
          !DO idx_comp=1,NUMBER_OF_COMPONENTS
          !   DO idx_dim=1,NUMBER_OF_DIMENSIONS
          !      IF(idx_dim==NUMBER_OF_DIMENSIONS) THEN
          !         LINE=LIST_STR(idx_comp)
          !      ELSE
          !         pos=INDEX(LINE,"*")                     
          !         LINE=EXTRACT(LIST_STR(idx_comp),1, pos-1)
          !         LIST_STR(idx_comp)=REMOVE(LIST_STR(idx_comp),1,pos)
          !      ENDIF
          !      CALL FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE(INTERPOLATION_XI(idx_comp, idx_dim), LINE, ERR, ERROR, *999)
          !   ENDDO
          !ENDDO          
       ENDIF !MASTER_COMPUTATIONAL_NUMBER/=my_computational_node_number
       
       !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LIST_ELEMENTAL_NODES:",LIST_ELEMENTAL_NODES(1),ERR,ERROR,*999)
       CALL MPI_BCAST(LIST_ELEMENTAL_NODES,number_of_node,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                                             
       CALL MPI_BCAST(LIST_COMP_NODAL_INDEX,number_of_node*NUMBER_OF_COMPONENTS,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                                      
       CALL MPI_BCAST(SHAPE_INDEX,SHAPE_SIZE,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
       CALL MPI_BCAST(LIST_COMP_NODES,NUMBER_OF_COMPONENTS,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
       CALL MPI_BCAST(MESH_COMPONENTS_OF_FIELD_COMPONENTS,NUMBER_OF_COMPONENTS,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
       CALL MPI_BCAST(INTERPOLATION_XI,NUMBER_OF_COMPONENTS*NUMBER_OF_DIMENSIONS,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,&
            &MPI_COMM_WORLD,MPI_IERROR)
       CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)                                             
       !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LIST_ELEMENTAL_NODES:",LIST_ELEMENTAL_NODES(1),ERR,ERROR,*999) 
       current_mesh_comp=1
       DO idx_comp=1, NUMBER_OF_COMPONENTS
          !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LIST_ELEMENTAL_NODES:",LIST_ELEMENTAL_NODES(1),ERR,ERROR,*999)
          IF(NUMBER_OF_DIMENSIONS==3) THEN
             CALL LIST_SEARCH(LIST_ELEMENT_NUMBER, SHAPE_INDEX(1),GLOBAL_ELEMENT_NUMBER, ERR,ERROR,*999)
          ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
             CALL LIST_SEARCH(LIST_ELEMENT_NUMBER, SHAPE_INDEX(2),GLOBAL_ELEMENT_NUMBER, ERR,ERROR,*999)
          ELSE IF(NUMBER_OF_DIMENSIONS==1) THEN
             CALL LIST_SEARCH(LIST_ELEMENT_NUMBER, SHAPE_INDEX(3),GLOBAL_ELEMENT_NUMBER, ERR,ERROR,*999)
          ELSE
             CALL FLAG_ERROR("Non recognized dimension size during reading elemental numbering",ERR,ERROR,*999)        
          ENDIF     
                    
          IF(MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp)==current_mesh_comp) THEN
             !find out whether the basis has been created
             pos=0
             DO idx_basis=1, BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS
                IF(SUM(BASIS_FUNCTIONS%BASES(idx_basis)%PTR%INTERPOLATION_XI(:)-INTERPOLATION_XI(idx_comp,:))==0) THEN
                   pos=idx_basis
                   EXIT
                ENDIF
             ENDDO
             
             IF(pos==0) THEN
                IF(ASSOCIATED(BASIS)) NULLIFY(BASIS)  
                CALL BASIS_CREATE_START(BASIS_FUNCTIONS%NUMBER_BASIS_FUNCTIONS+1,BASIS,ERR,ERROR,*999)
                CALL BASIS_NUMBER_OF_XI_SET(BASIS,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL BASIS_INTERPOLATION_XI_SET(BASIS,INTERPOLATION_XI(idx_comp,:),ERR,ERROR,*999)
                CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
             ENDIF
             
             CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(GLOBAL_ELEMENT_NUMBER,ELEMENTS_PTR(MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp))%PTR,&
                &BASIS,ERR,ERROR,*999)
             !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LIST_ELEMENTAL_NODES:",LIST_ELEMENTAL_NODES(1),ERR,ERROR,*999)                
             CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(GLOBAL_ELEMENT_NUMBER,ELEMENTS_PTR(MESH_COMPONENTS_OF_FIELD_COMPONENTS(idx_comp))%PTR,&
                &LIST_ELEMENTAL_NODES(LIST_COMP_NODAL_INDEX(idx_comp,:)),ERR,ERROR,*999)
             !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LIST_ELEMENTAL_NODES:",LIST_ELEMENTAL_NODES(1),ERR,ERROR,*999)   
             current_mesh_comp=current_mesh_comp+1
          ENDIF                   
       ENDDO !idx_comp    
       !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"FILE_END:",FILE_END,ERR,ERROR,*999)
       
    ENDDO !idx_exelems<NUMBER_OF_EXELEM_FILES
    
    DO idx_comp=1, NUMBER_OF_MESH_COMPONENTS
       CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,idx_comp, ERR,ERROR,*999)
    ENDDO           
    CALL MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*999)             
    
    IF(ALLOCATED(LIST_ELEMENT_NUMBER)) DEALLOCATE(LIST_ELEMENT_NUMBER)
    IF(ALLOCATED(ELEMENTS_PTR)) DEALLOCATE(ELEMENTS_PTR)
    !IF(ALLOCATED(LIST_NODAL_NUMBER)) DEALLOCATE(LIST_NODAL_NUMBER)
    IF(ALLOCATED(LIST_ELEMENTAL_NODES)) DEALLOCATE(LIST_ELEMENTAL_NODES)
    IF(ALLOCATED(LIST_STR)) DEALLOCATE(LIST_STR)
    IF(ALLOCATED(INTERPOLATION_XI)) DEALLOCATE(INTERPOLATION_XI)
    !IF(ALLOCATED(LIST_FIELD_COMPONENTS)) DEALLOCATE(LIST_FIELD_COMPONENTS)
    IF(ALLOCATED(MESH_COMPONENT_LOOKUP)) DEALLOCATE(MESH_COMPONENT_LOOKUP)
    IF(ALLOCATED(LIST_COMP_NODAL_INDEX)) DEALLOCATE(LIST_COMP_NODAL_INDEX)
    IF(ALLOCATED(LIST_COMP_NODES)) DEALLOCATE(LIST_COMP_NODES)
    IF(ASSOCIATED(BASIS)) NULLIFY(BASIS)       
                         
    CALL EXITS("FIELD_IO_IMPORT_GLOBAL_MESH")
    RETURN
999 CALL ERRORS("FIELD_IO_IMPORT_GLOBAL_MESH",ERR,ERROR)
    CALL EXITS("FIELD_IO_IMPORT_GLOBAL_MESH")
  END SUBROUTINE FIELD_IO_IMPORT_GLOBAL_MESH    

  !
  !================================================================================================================================
  !  
  
  !>Finding basis information 
  SUBROUTINE FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE(INTERPOLATION, LABEL_TYPE, ERR, ERROR, *)
    !Argument variables   
    INTEGER(INTG), INTENT(INOUT) :: INTERPOLATION !< xi interpolation type
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL_TYPE !<label type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   
    CALL ENTERS("FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE",ERR,ERROR,*999)    
       
    SELECT CASE(CHAR(LABEL_TYPE))
         CASE("l.Lagrange")
             INTERPOLATION=BASIS_LINEAR_LAGRANGE_INTERPOLATION
        CASE("q.Lagrange")
             INTERPOLATION=BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
         CASE("c.Lagrange")
             INTERPOLATION=BASIS_CUBIC_LAGRANGE_INTERPOLATION
         CASE("c.Hermite")
             INTERPOLATION=BASIS_CUBIC_HERMITE_INTERPOLATION
         CASE("q1.Hermite")
             INTERPOLATION=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
         CASE("q2.Hermite")
             INTERPOLATION=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
         CASE DEFAULT 
             CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
    END SELECT
                     
    CALL EXITS("FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE")
    RETURN
999 CALL ERRORS("FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE",ERR,ERROR)
    CALL EXITS("FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE")
  END SUBROUTINE FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE    

  !!
  !!================================================================================================================================
  !!
  !!>Create basis by reading information from multiple files in master node and broadcasting to others nodes    
  !SUBROUTINE FIELD_IO_CREATE_BASES_IN_LOCAL_PROCESS(NAME, BASES, MASTER_COMPUTATIONAL_NUMBER, & 
  !   & my_computational_node_number, ERR,ERROR,*)
  !  !Argument variables       
  !  TYPE(BASIS_FUNCTIONS_TYPE), POINTER :: BASES !<list of bases
  !  TYPE(VARYING_STRING), POINTER :: NAME !<name of input    
  !  INTEGER(INTG), INTENT(INOUT) :: MASTER_COMPUTATIONAL_NUMBER !< master number 
  !  INTEGER(INTG), INTENT(INOUT) :: my_computational_node_number !< my computational node number
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(BASIS_TYPE), POINTER :: BASIS
  !	TYPE(VARYING_STRING) :: LINE, CMISS_KEYWORD, FILE_NAME, FILE_STATUS
  !	TYPE(VARYING_STRING), ALLOCATABLE :: LIST_STR(:), LIST_STR1(:)
  !  INTEGER(INTG) :: FILE_ID !<file handle
  !  INTEGER(INTG) :: MPI_IERROR
  !  INTEGER(INTG) :: idx_basis, pos, NUM_INTERPOLATION_XI, NUMBER_OF_FILES, NUMBER_OF_BASES
  !  INTEGER(INTG) :: INTERPOLATION_XI(12), num_interp
  !  LOGICAL :: SWITCH_BASIS, FILE_EXIST
  !  
  !  CALL ENTERS("FIELD_IO_CREATE_BASES_IN_LOCAL_PROCESS", ERR,ERROR,*999)        
  !	                    
  !  !Intialise the bases
  !  idx_basis=0
  !  NUMBER_OF_BASES=0;
  !  ERR=0;
  !	
  !	IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN
  !    !the file name has to start from zero in a ascended order without break
  !    NUMBER_OF_FILES=0
  !    FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_FILES,"*",ERR,ERROR))//".exelem"
  !    INQUIRE(FILE=CHAR(FILE_NAME), EXIST=FILE_EXIST)
  !    IF(FILE_EXIST==.FALSE.) THEN
  !       CALL FLAG_ERROR("no file can be found, pls check again",ERR,ERROR,*999)
  !       GOTO 999                       
  !    ENDIF 
  !    DO WHILE(FILE_EXIST==.TRUE.)
  !       FILE_STATUS="OLD"
  !       FILE_ID=1030+NUMBER_OF_FILES
  !       CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)         
  !       CMISS_KEYWORD=", #Scale factor="
  !       DO WHILE(ERR==0)    
  !          CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
  !          IF(VERIFY(CMISS_KEYWORD,LINE)==0) THEN
  !             pos=INDEX(LINE,CMISS_KEYWORD)
  !             LINE=REMOVE(LINE, pos, LEN(LINE))
  !             IF(NUMBER_OF_BASES==0) THEN
  !                ALLOCATE(LIST_STR(1), STAT=ERR)
  !                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate character buffer in reading basis",ERR,ERROR,*999)
  !                LIST_STR(1)=LINE
  !             ELSE
  !                SWITCH_BASIS=.FALSE.         
  !                DO idx_basis=1,NUMBER_OF_BASES
  !                   IF(LIST_STR(idx_basis)==LINE) THEN
  !                      SWITCH_BASIS=.TRUE.
  !                      EXIT
  !                   ENDIF        
  !                ENDDO 
  !                IF(SWITCH_BASIS==.FALSE.) THEN
  !                   ALLOCATE(LIST_STR1(NUMBER_OF_BASES), STAT=ERR)
  !                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate character buffer in reading basis",ERR,ERROR,*999)
  !                   LIST_STR1(:)=LIST_STR(:)
  !                   DEALLOCATE(LIST_STR)
  !                   ALLOCATE(LIST_STR(NUMBER_OF_BASES+1), STAT=ERR)
  !                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate character buffer in reading basis",ERR,ERROR,*999)
  !                   LIST_STR(1:NUMBER_OF_BASES)=LIST_STR1(:)
  !                   LIST_STR(NUMBER_OF_BASES+1)=LINE
  !                   NUMBER_OF_BASES=NUMBER_OF_BASES+1   
  !                   DEALLOCATE(LIST_STR1)             
  !                ENDIF !SWITCH_BASIS==.FALSE.
  !             ENDIF !NUMBER_OF_BASES==0     
  !          ENDIF  !VERIFY(CMISS_KEYWORD,LINE)==0 
  !       ENDDO !(ERR==0)                 
  !       
  !       CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
  !       !checking the next file
  !       NUMBER_OF_FILES=NUMBER_OF_FILES+1
  !       FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_FILES,"*",ERR,ERROR))//".exelem"
  !       INQUIRE(FILE=CHAR(FILE_NAME), EXIST=FILE_EXIST)
  !    ENDDO!FILE_EXIST==.TRUE.
  !    !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Total number of exelment files = ",NUMBER_OF_FILES-1, ERR,ERROR,*999)
  !  ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number
  !  
  !  !broadcasting the number of bases
  !  CALL MPI_BCAST(NUMBER_OF_BASES,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
  !  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
  !  
  !  IF(NUMBER_OF_BASES/=0) THEN
  !     ALLOCATE(BASES%BASES(NUMBER_OF_BASES), STAT=ERR)
  !     IF(ERR/=0) CALL FLAG_ERROR("Could not allocate basis buffer in reading basis",ERR,ERROR,*999)
  !     BASES%NUMBER_BASIS_FUNCTIONS=NUMBER_OF_BASES
  !  ELSE
  !     CALL FLAG_ERROR("can not file any basis informations in all exelem files",ERR,ERROR,*999)
  !     GOTO 999                                 
  !  ENDIF
  !  !BASES%BASES=>BASIS_PTR
  !  NUM_INTERPOLATION_XI=0
  !  DO idx_basis=1,NUMBER_OF_BASES
  !     CALL BASIS_CREATE_START(idx_basis,BASIS,ERR,ERROR,*999)
  !     
  !     IF(MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number) THEN
  !        DO WHILE(VERIFY("*",LIST_STR(idx_basis))==0)
  !           NUM_INTERPOLATION_XI=NUM_INTERPOLATION_XI+1 
  !           pos=INDEX(LIST_STR(idx_basis),"*")
  !           LINE=EXTRACT(LIST_STR(idx_basis), 1, pos)
  !           LIST_STR(idx_basis)=REMOVE(LIST_STR(idx_basis),1,pos)
  !           CALL FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE(num_interp, LINE, ERR, ERROR, *999)
  !           INTERPOLATION_XI(NUM_INTERPOLATION_XI)=num_interp
  !        ENDDO
  !        NUM_INTERPOLATION_XI=NUM_INTERPOLATION_XI+1
  !        LINE=LIST_STR(idx_basis)
  !        CALL FIELD_IO_TRANSLATE_LABEL_INTO_INTERPOLATION_TYPE(num_interp, LINE, ERR, ERROR, *999)
  !        INTERPOLATION_XI(NUM_INTERPOLATION_XI)=num_interp
  !     ENDIF !MASTER_COMPUTATIONAL_NUMBER==my_computational_node_number       
  !
  !     !broadcasting the number of XI
  !     CALL MPI_BCAST(NUM_INTERPOLATION_XI,1,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
  !     CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
  !     CALL BASIS_NUMBER_OF_XI_SET(BASIS,NUM_INTERPOLATION_XI,ERR,ERROR,*999)
  !
  !     !broadcasting the number of XI
  !     CALL MPI_BCAST(INTERPOLATION_XI,NUM_INTERPOLATION_XI,MPI_INTEGER,MASTER_COMPUTATIONAL_NUMBER,MPI_COMM_WORLD,MPI_IERROR)
  !     CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)       
  !     CALL BASIS_INTERPOLATION_XI_SET(BASIS,INTERPOLATION_XI,ERR,ERROR,*999)
  !
  !     CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
  !     
  !     BASES%BASES(idx_basis)%PTR=>BASIS      
  !  ENDDO !idx_basis
  ! 
  !  IF(ALLOCATED(LIST_STR)) DEALLOCATE(LIST_STR)
  !  IF(ALLOCATED(LIST_STR1)) DEALLOCATE(LIST_STR1)
  !  
  !  CALL EXITS("FIELD_IO_CREATE_BASES_IN_LOCAL_PROCESS")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_CREATE_BASES_IN_LOCAL_PROCESS",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_CREATE_BASES_IN_LOCAL_PROCESS")
  !  RETURN 1  
  !END SUBROUTINE FIELD_IO_CREATE_BASES_IN_LOCAL_PROCESS


  !
  !================================================================================================================================
  !
  
  !>Export elemental information into multiple files \see{FIELD_IO::FIELD_IO_ELEMENTS_EXPORT}.     
  SUBROUTINE FIELD_IO_ELEMENTS_EXPORT(FIELDS, FILE_NAME, METHOD,ERR,ERROR,*)
  !checking the input data for IO and initialize the nodal information set
  !the following items will be checked: the region (the same?), all the pointer(valid?)   
  !in this version, different decomposition method will be allowed for the list of field variables(but still in the same region)
  !even the each process has exactly the same nodal information and each process will write out exactly the same data
  !because CMGui can read the same data for several times   

    !Argument variables       
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<the field object
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<file name
    TYPE(VARYING_STRING), INTENT(IN):: METHOD !< method used for IO
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET) :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<elemental information in this process
    INTEGER(INTG):: my_computational_node_number !<local process number    
    INTEGER(INTG):: computational_node_numbers   !<total process numbers      

    CALL ENTERS("FIELD_IO_ELEMENTS_EXPORT", ERR,ERROR,*999)    
       
    !Get the number of computational nodes
    computational_node_numbers=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !Get my computational node number
    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    IF(METHOD=="FORTRAN") THEN
       CALL FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
       CALL FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*999)
       CALL FIELD_IO_ELEMENTAL_INFO_SET_SORT(LOCAL_PROCESS_ELEMENTAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)    
       CALL FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FILE_NAME, my_computational_node_number, &
            &computational_node_numbers, ERR, ERROR, *999)
       CALL FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*999)
    ELSE IF(METHOD=="MPIIO") THEN
       CALL FLAG_ERROR("MPI IO has not been implemented yet",ERR,ERROR,*999)
    ELSE 
       CALL FLAG_ERROR("Unknown method!",ERR,ERROR,*999)   
    ENDIF   
     
    CALL EXITS("FIELD_IO_ELEMENTS_EXPORT")
    RETURN
999 CALL ERRORS("FIELD_IO_ELEMENTS_EXPORT",ERR,ERROR)
    CALL EXITS("FIELD_IO_ELEMENTS_EXPORT")
    RETURN 1  
  END SUBROUTINE FIELD_IO_ELEMENTS_EXPORT
  
  !
  !================================================================================================================================
  !  
  
  !>Finding basis information 
  FUNCTION FIELD_IO_BASIS_LHTP_FAMILY_LABEL(BASIS, num_scl, num_node, LABEL_TYPE, ERR, ERROR)
    !Argument variables   
    TYPE(BASIS_TYPE), INTENT(IN) :: BASIS !<The error string
    INTEGER(INTG), INTENT(INOUT) :: num_scl
    INTEGER(INTG), INTENT(INOUT) :: num_node
    INTEGER(INTG), INTENT(IN) ::LABEL_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) ::FIELD_IO_BASIS_LHTP_FAMILY_LABEL
    INTEGER(INTG) :: ni 
   
    CALL ENTERS("FIELD_IO_BASIS_LHTP_FAMILY_LABEL",ERR,ERROR,*999)    
   
    IF(BASIS%NUMBER_OF_XI==0) CALL FLAG_ERROR("number of xi in the basis is zero",ERR,ERROR,*999)           
    
    num_scl=1;
    num_node=1
    DO ni=1,BASIS%NUMBER_OF_XI-1
       SELECT CASE(BASIS%INTERPOLATION_XI(ni))
         CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"l.Lagrange*"
             num_scl=num_scl*2  
             num_node=num_node*2
          CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"q.Lagrange*"
             num_scl=num_scl*3
             num_node=num_node*3
          CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"c.Lagrange*"
             num_scl=num_scl*4
             num_node=num_node*4
          CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"c.Hermite*"
             num_scl=num_scl*2*2
             num_node=num_node*2
          CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"q1.Hermite*"
             num_scl=num_scl*2*2
             num_node=num_node*2            
          CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"q2.Hermite*"
             num_scl=num_scl*2*2
             num_node=num_node*2
          CASE DEFAULT 
             CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
       END SELECT
       !IF(BASIS%COLLAPSED_XI(ni)==BASIS_XI_COLLAPSED) THEN
       !   BASIS%NUMBER_OF_COLLAPSED_XI=BASIS%NUMBER_OF_COLLAPSED_XI+1
       !   COLLAPSED_XI(BASIS%NUMBER_OF_COLLAPSED_XI)=ni
       !   BASIS%DEGENERATE=.TRUE.
       !ENDIF
       !NUMBER_OF_NODES=NUMBER_OF_NODES*BASIS%NUMBER_OF_NODES_XI(ni)
       !IF(BASIS%NUMBER_OF_NODES_XI(ni)>MAX_NUM_NODES) MAX_NUM_NODES=BASIS%NUMBER_OF_NODES_XI(ni)
    ENDDO !ni
    DO ni=BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
       SELECT CASE(BASIS%INTERPOLATION_XI(ni))
         CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"l.Lagrange"
             num_scl=num_scl*2  
             num_node=num_node*2
          CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"q.Lagrange"
             num_scl=num_scl*3
             num_node=num_node*3
          CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"c.Lagrange"
             num_scl=num_scl*4
             num_node=num_node*4
          CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"c.Hermite"
             num_scl=num_scl*2*2
             num_node=num_node*2
          CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"q1.Hermite"
             num_scl=num_scl*2*2
             num_node=num_node*2            
          CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
             FIELD_IO_BASIS_LHTP_FAMILY_LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL//"q2.Hermite"
             num_scl=num_scl*2*2
             num_node=num_node*2
          CASE DEFAULT 
             CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
       END SELECT
       !IF(BASIS%COLLAPSED_XI(ni)==BASIS_XI_COLLAPSED) THEN
       !   BASIS%NUMBER_OF_COLLAPSED_XI=BASIS%NUMBER_OF_COLLAPSED_XI+1
       !   COLLAPSED_XI(BASIS%NUMBER_OF_COLLAPSED_XI)=ni
       !   BASIS%DEGENERATE=.TRUE.
       !ENDIF
       !NUMBER_OF_NODES=NUMBER_OF_NODES*BASIS%NUMBER_OF_NODES_XI(ni)
       !IF(BASIS%NUMBER_OF_NODES_XI(ni)>MAX_NUM_NODES) MAX_NUM_NODES=BASIS%NUMBER_OF_NODES_XI(ni)
    ENDDO !ni

    !FIELD_IO_BASIS_LHTP_FAMILY_LABEL=" "!FIELD_IO_BASIS_LHTP_FAMILY_LABEL(1:(LEN_TRIM(FIELD_IO_BASIS_LHTP_FAMILY_LABEL)-1))
    SELECT CASE(LABEL_TYPE)
       CASE (FIELD_IO_SCALE_FACTORS_NUMBER_TYPE)
          FIELD_IO_BASIS_LHTP_FAMILY_LABEL=TRIM(FIELD_IO_BASIS_LHTP_FAMILY_LABEL)//", #Scale factors="//TRIM(NUMBER_TO_VSTRING&
               &(num_scl,"*",ERR,ERROR))
       CASE (FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE)
          FIELD_IO_BASIS_LHTP_FAMILY_LABEL=TRIM(FIELD_IO_BASIS_LHTP_FAMILY_LABEL)//", no modify, standard node based"
       CASE DEFAULT
          CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
    END SELECT   
    
    !MAX_SCALE_FACTORS=MAX(MAX_SCALE_FACTORS,num_scl)
                     
    CALL EXITS("FIELD_IO_BASIS_LHTP_FAMILY_LABEL")
    RETURN
999 CALL ERRORS("FIELD_IO_BASIS_LHTP_FAMILY_LABEL",ERR,ERROR)
    CALL EXITS("FIELD_IO_BASIS_LHTP_FAMILY_LABEL")
  END FUNCTION FIELD_IO_BASIS_LHTP_FAMILY_LABEL    

  !
  !================================================================================================================================
  !  

  !>Write the header of a group elements using FORTRAN 
  SUBROUTINE FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN(LOCAL_PROCESS_ELEMENTAL_INFO_SET, LOCAL_ELEMENTAL_NUMBER,  MAX_NODE_CPMP_INDEX,&
          &NUM_OF_SCALING_FACTOR_SETS, LIST_COMP_SCALE, my_computational_node_number, FILE_ID, ERR,ERROR, *)
    !Argument variables  
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_ELEMENTAL_INFO_SET  !<LOCAL_PROCESS_NODAL_INFO_SET           
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENTAL_NUMBER !<element number
    INTEGER(INTG), INTENT(INOUT) ::  MAX_NODE_CPMP_INDEX !<MAX_NODE_INDEX
    INTEGER(INTG), INTENT(INOUT) :: NUM_OF_SCALING_FACTOR_SETS !<NUM_OF_SCALING_FACTOR_SETS
    INTEGER(INTG), INTENT(INOUT) :: LIST_COMP_SCALE(:)
    INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN) :: FILE_ID !< FILE ID    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: field_ptr
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable_ptr   
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_ELEMENTS !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS ! domain nodes
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
    TYPE(BASIS_TYPE), POINTER :: BASIS 
    TYPE(BASIS_TYPE), POINTER :: BASIS1
    TYPE(VARYING_STRING) :: LINE, LABEL
    INTEGER(INTG), ALLOCATABLE :: LIST_SCALE_FACTORS(:), GROUP_LOCAL_NUMBER(:), GROUP_SCALE_FACTORS(:)
    INTEGER(INTG), ALLOCATABLE :: GROUP_NODE(:), GROUP_VARIABLES(:)
    INTEGER(INTG) :: NUM_OF_VARIABLES, MAX_NUM_NODES !NUM_OF_FIELDS NUM_OF_NODES
    INTEGER(INTG) :: domain_idx, domain_no, MY_DOMAIN_INDEX, local_number, global_number
    INTEGER(INTG) ::nn, np, mk, nk, num_scl, num_node, comp_idx, scl_idx, scl_idx1, var_idx!value_idx field_idx global_var_idx comp_idx1 ny2 
    LOGICAL :: SWITCH   
   
    CALL ENTERS("FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN",ERR,ERROR,*999)    
    
    !colllect nodal header information for IO first
    
    !!get the number of this computational node from mpi pool
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
    !IF(ERR/=0) GOTO 999        
    
    !attach the temporary pointer
    !tmp_components=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS    
    
    !collect maximum number of nodal derivatives, number of fields and variables 
    NUM_OF_SCALING_FACTOR_SETS=0
    NUM_OF_VARIABLES=0
    MAX_NUM_NODES=0
    MAX_NODE_CPMP_INDEX=0
    global_number=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_ELEMENTAL_NUMBER)
    NULLIFY(field_ptr)     
    NULLIFY(variable_ptr)
    IF(ALLOCATED(GROUP_LOCAL_NUMBER)) DEALLOCATE(GROUP_LOCAL_NUMBER)
    ALLOCATE(GROUP_LOCAL_NUMBER(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS),&
        &STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate GROUP_LOCAL_NUMBER in exelem header",ERR,ERROR,*999)
    IF(ALLOCATED(LIST_SCALE_FACTORS)) DEALLOCATE(LIST_SCALE_FACTORS)
    ALLOCATE(LIST_SCALE_FACTORS(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS),&
        &STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate LIST_SCALE_FACTORS in exelem header",ERR,ERROR,*999)
     
    DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of fields
       !IF (.NOT.ASSOCIATED(field_ptr, target=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD)) THEN
       !      NUM_OF_FIELDS=NUM_OF_FIELDS+1
       !    field_ptr=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD
       !ENDIF         
       
       !calculate the number of variables
       IF (.NOT.ASSOCIATED(variable_ptr, target=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE)) THEN
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          variable_ptr=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%&
              &PTR%FIELD_VARIABLE
       ENDIF      
              
       !finding the local numbering through the global to local mapping   
       DOMAIN_MAPPING_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%&
           &PTR%DOMAIN%MAPPINGS%ELEMENTS       
       !get the domain index for this variable component according to my own computional node number
       DO domain_idx=1,DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          IF(domain_no==my_computational_node_number) THEN
             MY_DOMAIN_INDEX=domain_idx
             EXIT !out of loop--domain_idx
          ENDIF
       ENDDO !domain_idX
       local_number=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
       GROUP_LOCAL_NUMBER(comp_idx)=local_number
       !use local domain information find the out the maximum number of derivatives
       DOMAIN_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%&
           &DOMAIN%TOPOLOGY%ELEMENTS 
       DOMAIN_NODES=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%&
           &DOMAIN%TOPOLOGY%NODES
       BASIS=>DOMAIN_ELEMENTS%ELEMENTS(local_number)%BASIS
       IF(BASIS%NUMBER_OF_NODES>MAX_NUM_NODES) THEN
          MAX_NODE_CPMP_INDEX=comp_idx
          MAX_NUM_NODES=BASIS%NUMBER_OF_NODES
       ENDIF   
       IF(BASIS%DEGENERATE==.FALSE.)  THEN
          IF(comp_idx==1) THEN
             NUM_OF_SCALING_FACTOR_SETS=NUM_OF_SCALING_FACTOR_SETS+1
             LIST_SCALE_FACTORS(NUM_OF_SCALING_FACTOR_SETS)=comp_idx
          ELSE
             SWITCH=.FALSE.
             DO scl_idx1=1, NUM_OF_SCALING_FACTOR_SETS
                BASIS1=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%&
                   &COMPONENTS(LIST_SCALE_FACTORS(scl_idx1))%PTR%DOMAIN%TOPOLOGY%ELEMENTS% &
                & ELEMENTS(GROUP_LOCAL_NUMBER(LIST_SCALE_FACTORS(scl_idx1)))%BASIS
                IF(SUM(BASIS1%INTERPOLATION_XI(1:BASIS1%NUMBER_OF_XI)-BASIS%INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI))==0.AND.&
                   &SUM(BASIS1%INTERPOLATION_TYPE(1:BASIS1%NUMBER_OF_XI)-BASIS%INTERPOLATION_TYPE(1:BASIS%NUMBER_OF_XI))==0.AND. &
                   &SUM(BASIS1%INTERPOLATION_ORDER(1:BASIS1%NUMBER_OF_XI)-BASIS%INTERPOLATION_ORDER(1:BASIS%NUMBER_OF_XI))==0) THEN
                   SWITCH=.TRUE.
                   EXIT
                ENDIF                   
             ENDDO!scl_idx1
             IF(SWITCH==.FALSE.) THEN
                NUM_OF_SCALING_FACTOR_SETS=NUM_OF_SCALING_FACTOR_SETS+1
                LIST_SCALE_FACTORS(NUM_OF_SCALING_FACTOR_SETS)=comp_idx
             ENDIF                 
          ENDIF
          LIST_COMP_SCALE(comp_idx)=NUM_OF_SCALING_FACTOR_SETS
       ENDIF!BASIS%DEGENERATE=.FALSE.
    ENDDO !comp_idx         
    !!Allocate the momery for group of field variables
    !ALLOCATE(GROUP_FIELDS(NUM_OF_FIELDS),STAT=ERR)
    !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty field buffer in IO",ERR,ERROR,*999)
    !!Allocate the momery for group of field components   
    IF(ALLOCATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)    
    ALLOCATE(GROUP_VARIABLES(NUM_OF_VARIABLES),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
    
    !!Allocate the momery for group of maximum number of derivatives     
    IF(ALLOCATED(GROUP_SCALE_FACTORS)) DEALLOCATE(GROUP_SCALE_FACTORS)    
    ALLOCATE(GROUP_SCALE_FACTORS(NUM_OF_SCALING_FACTOR_SETS),STAT=ERR)    
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
    
    IF(ALLOCATED(GROUP_NODE)) DEALLOCATE(GROUP_NODE)    
    ALLOCATE(GROUP_NODE(NUM_OF_SCALING_FACTOR_SETS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)            
        
    !fill information into the group of fields and variables
    NULLIFY(variable_ptr) 
    GROUP_VARIABLES(:)=0         
    NUM_OF_VARIABLES=0
    DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of variables
       IF (.NOT.ASSOCIATED(variable_ptr, target=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%&
          &COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE)) THEN
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          variable_ptr=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%&
             &FIELD_VARIABLE
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
       ELSE
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1   
       ENDIF                    
    ENDDO  !comp_idx                
    
    !write out the scale factor set information
    LINE=" #Scale factor sets= "//TRIM(NUMBER_TO_VSTRING(NUM_OF_SCALING_FACTOR_SETS,"*",ERR,ERROR))
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    !MAX_SCALE_FACTORS=0
    GROUP_SCALE_FACTORS(:)=0
    DO scl_idx=1,NUM_OF_SCALING_FACTOR_SETS                       
       BASIS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(LIST_SCALE_FACTORS(scl_idx))% &
         &PTR%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(GROUP_LOCAL_NUMBER(LIST_SCALE_FACTORS(scl_idx)))%BASIS
       num_scl=0
       IF(ASSOCIATED(BASIS)) THEN
          SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
             LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL(BASIS, num_scl, num_node, FIELD_IO_SCALE_FACTORS_NUMBER_TYPE, ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get basis type of lagrange_hermite label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF               
            !CASE(BASIS_SIMPLEX_TYPE)
            !  CALL BASIS_SIMPLEX_FAMILY_CREATE(BASIS,ERR,ERROR,*999)
            CASE DEFAULT
               CALL FLAG_ERROR("Basis type "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))//" is invalid or not implemented",&
               &ERR,ERROR,*999)
          END SELECT
       ELSE
          CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
       ENDIF
       GROUP_SCALE_FACTORS(scl_idx)=num_scl!numer of scale factors in scale factor set
       GROUP_NODE(scl_idx)=num_node!numer of nodes in scale factor set
       LINE="   "//LABEL
       CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                                    
    ENDDO!scl_idx

    LINE=" #Nodes=           "//TRIM(NUMBER_TO_VSTRING(MAX_NUM_NODES,"*",ERR,ERROR))
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)      
    LINE=" #Fields="//TRIM(NUMBER_TO_VSTRING(NUM_OF_VARIABLES,"*",ERR,ERROR))
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)      
        
    !write out the nodal header
    var_idx=0
    NULLIFY(variable_ptr) 
    !comp_idx=1
    !field_idx=1
    !value_idx=1
    !comp_idx1=1
    !global_var_idx=0
    !NUM_OF_VARIABLES=0
    DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of fields
       !IF (.NOT.ASSOCIATED(field_ptr, target=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD)) THEN
       !      NUM_OF_FIELDS=NUM_OF_FIELDS+1
       !    field_ptr=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD
       !ENDIF         
       
       !calculate the number of variables
       !grouping field variables and components together
       IF(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%&
          & COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE)) THEN !different variables            
          var_idx=var_idx+1
          variable_ptr=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%&
               &PTR%FIELD_VARIABLE          
          !write out the field information
          LABEL=" "//TRIM(NUMBER_TO_VSTRING(var_idx,"*",ERR,ERROR))//") "&
          &//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%&
              &COMPONENTS(comp_idx)%PTR, FIELD_IO_VARIABLE_LABEL,ERR,ERROR)
          IF(ERR/=0) THEN
             CALL FLAG_ERROR("can not get variable label",ERR,ERROR,*999)     
             GOTO 999               
          ENDIF        
          LINE=TRIM(LABEL)//", #Components="//TRIM(NUMBER_TO_VSTRING(GROUP_VARIABLES(var_idx),"*",ERR,ERROR))                  
          CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
       ENDIF

       !write out the component information
       LABEL="   "//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%&
            &COMPONENTS(comp_idx)%PTR, FIELD_IO_COMPONENT_LABEL,ERR,ERROR)
       IF(ERR/=0) THEN
          CALL FLAG_ERROR("can not get component label",ERR,ERROR,*999)     
          GOTO 999               
       ENDIF        
       LINE=TRIM(LABEL)//"."                          
       BASIS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%&
           &TOPOLOGY%ELEMENTS%ELEMENTS(comp_idx)%BASIS
       SELECT CASE(BASIS%TYPE)
         CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            LABEL=FIELD_IO_BASIS_LHTP_FAMILY_LABEL(BASIS, num_scl, num_node, FIELD_IO_SCALE_FACTORS_PROPERTY_TYPE, ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get basis type of lagrange_hermite label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF                           
            !CASE(BASIS_SIMPLEX_TYPE)
            !  CALL BASIS_SIMPLEX_FAMILY_CREATE(BASIS,ERR,ERROR,*999)
         CASE DEFAULT
            CALL FLAG_ERROR("Basis type "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))//" is invalid or not implemented",&
                 & ERR,ERROR,*999)
       END SELECT
       LINE=TRIM(LINE)//"  "//TRIM(LABEL) 
       CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999) 
       
       LINE="     #Nodes= "//TRIM(NUMBER_TO_VSTRING(num_node,"*",ERR,ERROR))
       CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)      

       IF(BASIS%DEGENERATE==.FALSE.) THEN
          IF(LIST_COMP_SCALE(comp_idx)==1) THEN
             scl_idx=0
          ELSE   
             scl_idx= SUM(GROUP_SCALE_FACTORS(1:LIST_COMP_SCALE(comp_idx)-1))
          ENDIF   
          BASIS=>DOMAIN_ELEMENTS%ELEMENTS(LOCAL_ELEMENTAL_NUMBER)%BASIS    
          DO nn=1,BASIS%NUMBER_OF_NODES
             LINE="      "//TRIM(NUMBER_TO_VSTRING(nn,"*",ERR,ERROR))//".  #Values="//&
                 &TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_DERIVATIVES(nn),"*",ERR,ERROR))
             CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                     
             np=DOMAIN_ELEMENTS%ELEMENTS(LOCAL_ELEMENTAL_NUMBER)%ELEMENT_NODES(nn)               
             LINE="       Value indices:  "             
             DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                nk=DOMAIN_ELEMENTS%ELEMENTS(LOCAL_ELEMENTAL_NUMBER)%ELEMENT_DERIVATIVES(mk,nn)
                !ny2=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk)
                !GROUP_DERIVATIVES(mk)=ny2
                LINE=LINE//TRIM(NUMBER_TO_VSTRING(nk, "*",ERR,ERROR))
             ENDDO !mk
             CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                 
  
             LINE="       Scale factor indices: "             
             DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                scl_idx=scl_idx+1
                LINE=LINE//TRIM(NUMBER_TO_VSTRING(scl_idx, "*",ERR,ERROR))
             ENDDO !mk
               CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                                         
          ENDDO !nn
       ELSE
          CALL FLAG_ERROR("exporting degenerated nodes has not been implemented",ERR,ERROR,*999) 
       ENDIF  
    ENDDO !comp_idx         
    
    
    !release temporary memory
    IF(ALLOCATED(LIST_SCALE_FACTORS)) DEALLOCATE(LIST_SCALE_FACTORS)
    IF(ALLOCATED(GROUP_LOCAL_NUMBER)) DEALLOCATE(GROUP_LOCAL_NUMBER)
    IF(ALLOCATED(GROUP_SCALE_FACTORS)) DEALLOCATE(GROUP_SCALE_FACTORS)    
    IF(ALLOCATED(GROUP_NODE)) DEALLOCATE(GROUP_NODE)    
    IF(ALLOCATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)
    !IF(ASSOCIATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)     
        
    CALL EXITS("FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN")
    RETURN
999 CALL ERRORS("FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN",ERR,ERROR)
    CALL EXITS("FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN")
    RETURN 1       
  END SUBROUTINE FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN  
  
  !
  !================================================================================================================================
  !  

  !>Write all the elemental information from LOCAL_PROCESS_NODAL_INFO_SET to exelem files 
  SUBROUTINE FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, NAME, my_computational_node_number, &
  &computational_node_numbers,ERR, ERROR, *)
    !the reason that my_computational_node_number is used in the argument is for future extension
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<nodal information in this process
    TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN):: computational_node_numbers !<total process number       
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: FILE_NAME, LINE, FILE_STATUS !the prefix name of file.    
    TYPE(BASIS_TYPE), POINTER ::BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_ELEMENTS !The domain mapping to calculate elemental mappings
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS ! domain elements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain elements    
    INTEGER(INTG) :: FILE_ID, domain_idx, domain_no, local_number, global_number, MY_DOMAIN_INDEX, MAX_NODE_CPMP_INDEX, NUM_DIM
    INTEGER(INTG), ALLOCATABLE :: LIST_COMP_SCALE(:), NODAL_NUMBER(:)!LIST_COMP(:) !Components which will be used for export scale factors    
    INTEGER(INTG) :: nk, np, nn, ns, mk, ny2, elem_idx, comp_idx,  scal_idx, NUM_OF_SCALING_FACTOR_SETS !dev_idx  elem_num  
    REAL(DP), ALLOCATABLE :: SCALE_FACTORS(:), SCALING_BUFFER(:) 

    CALL ENTERS("FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE",ERR,ERROR,*999)    
    
    !get my own computianal node number--be careful the rank of process in the MPI pool 
    !is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
    !will be a secured way to get the number         
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)         
    !IF(ERR/=0) GOTO 999
    FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(my_computational_node_number,"*",ERR,ERROR))//".exelem"        
    FILE_ID=1029+my_computational_node_number
    NUM_OF_SCALING_FACTOR_SETS=0
    
    IF(.NOT.ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)) THEN
       CALL FLAG_ERROR("the elemental information set in input is invalid",ERR,ERROR,*999)            
    ENDIF

    IF(.NOT.ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("the elemental information set is not associated with any numbering list",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS==0) THEN 
       CALL FLAG_ERROR("the elemental information set does not contain any nodes",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(1)%SAME_HEADER==.TRUE.) THEN
       CALL FLAG_ERROR("the first header flag of elemental information set should be false",ERR,ERROR,*999)            
    ENDIF
    
    !NULLIFY(SCALE_FACTORS)    
    !NULLIFY(SCALING_BUFFER)
    !NULLIFY(LIST_COMP_SCALE)
    !NULLIFY(tmp_components)
    
    NUM_DIM=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(1)%COMPONENTS(1)%PTR%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
    
    !open a file 
    FILE_STATUS="REPLACE"
    CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
    
    !write out the group name    
    LINE=FIELD_IO_FILEDS_GROUP_INFO_GET(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS, ERR,ERROR)    
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get group name in IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF               
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    !write out the number of files    
    LINE=FIELD_IO_MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)      
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get multiple file information in IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF             
    !CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_Thttp://www.bioeng.auckland.ac.nz/home/home.phpRIM(LINE), ERR,ERROR,*999)       

    !write out elements
    LINE=" Shape.  Dimension="//TRIM(NUMBER_TO_VSTRING(NUM_DIM,"*", ERR,ERROR))
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
     
    ns=1!
    DO elem_idx=1, LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
      
       !tmp_components=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS         
       global_number=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(elem_idx)
       

       IF(ALLOCATED(LIST_COMP_SCALE)) DEALLOCATE(LIST_COMP_SCALE)         
       ALLOCATE(LIST_COMP_SCALE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%NUMBER_OF_COMPONENTS),STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate LIST_COMP_SCALE in exelem io",ERR,ERROR,*999)

       !check whether need to write out the nodal information header  
       IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%SAME_HEADER==.FALSE.) THEN
          !print * ,"elem_idx:", elem_idx, "in", my_computational_node_number
           !write out the nodal header
           CALL FIELD_IO_EXPORT_ELEMENTAL_GROUP_HEADER_FORTRAN(LOCAL_PROCESS_ELEMENTAL_INFO_SET, elem_idx, MAX_NODE_CPMP_INDEX, &
              &NUM_OF_SCALING_FACTOR_SETS, &
              &LIST_COMP_SCALE ,my_computational_node_number, FILE_ID, ERR,ERROR,*999) 
           !value_idx=value_idx-1 !the len of NODAL_BUFFER                      
          !checking: whether need to allocate temporary memory for Io writing
       ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(elem_idx)%SAME_HEADER==.FALSE.

       !IF(ASSOCIATED(SCALING_BUFFER)) THEN
       !   IF(SIZE(SCALING_BUFFER)<MAX_NUM_OF_SCALING_INDICES) THEN 
       !      DEALLOCATE(SCALING_BUFFER)
       !      !DEALLOCATE(GROUP_DERIVATIVES)
       !      ALLOCATE(SCALING_BUFFER(MAX_NUM_OF_SCALING_INDICES),STAT=ERR)
       !      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty scaling buffer in IO writing",ERR,ERROR,*999)
       !      !ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
       !      !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
       !   ENDIF
       !ELSE
       !   ALLOCATE(SCALING_BUFFER(MAX_NUM_OF_SCALING_INDICES),STAT=ERR)
       !   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty scaling buffer in IO writing",ERR,ERROR,*999)
       !   !ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
       !   !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
       !ENDIF
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !write out elemental information        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !element info              
       IF(NUM_DIM==3) THEN
          LINE=" Element:            "//TRIM(NUMBER_TO_VSTRING(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%&
               &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)% &
               &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%USER_NUMBER,"*", &
               &ERR,ERROR))//" 0 0"
       ELSE IF(NUM_DIM==2) THEN
          LINE=" Element:            "//"0 "//TRIM(NUMBER_TO_VSTRING(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%&
               &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)% &
               &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%USER_NUMBER,"*", &
               &ERR,ERROR))//" 0"
       
       ELSE IF(NUM_DIM==1) THEN
          LINE=" Element:            "//" 0 0"//TRIM(NUMBER_TO_VSTRING(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%&
               &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)% &
               &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%USER_NUMBER,"*", &
               &ERR,ERROR))       
       ELSE
          CALL FLAG_ERROR("non-recognized dimension size in exporting exelem files",ERR,ERROR,*999)   
       ENDIF        
       CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
       !node info
       LINE="   Nodes:" !NEED TO LIST ALL THE NODES(SUCH CUBIC CASE)
       CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)
       BASIS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%&
          &TOPOLOGY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%&
          &MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%BASIS
       IF(BASIS%DEGENERATE==.FALSE.) THEN
          ALLOCATE(NODAL_NUMBER(BASIS%NUMBER_OF_NODES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list of nodal numbering when writing out exelem",ERR,ERROR,*999)                
          
          DO nn=1,BASIS%NUMBER_OF_NODES
             NODAL_NUMBER(nn)=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%&
                 &COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%&
                 &ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR&
                 &%ELEMENTS%ELEMENTS(global_number)%USER_ELEMENT_NODES(nn)    
          ENDDO !nn
          CALL FIELD_IO_FORTRAN_FILE_WRITE_INTG(FILE_ID, NODAL_NUMBER, BASIS%NUMBER_OF_NODES, ERR,ERROR,*999)
          DEALLOCATE(NODAL_NUMBER)
       ELSE
          CALL FLAG_ERROR("exporting degenerated nodes has not been implemented",ERR,ERROR,*999) 
       ENDIF       
                         
       !write out scale factors information        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       LINE="   Scale factors:"
       CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       scal_idx=1
       DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%NUMBER_OF_COMPONENTS          

          !finding the local numbering through the global to local mapping   
          DOMAIN_MAPPING_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(comp_idx)%PTR%&
             &DOMAIN%MAPPINGS%ELEMENTS       
          !get the domain index for this variable component according to my own computional node number
          DO domain_idx=1,DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
             domain_no=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
             IF(domain_no==my_computational_node_number) THEN
                MY_DOMAIN_INDEX=domain_idx
                EXIT !out of loop--domain_idx
             ENDIF
          ENDDO !domain_idX
          local_number=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
          !use local domain information find the out the maximum number of derivatives
          DOMAIN_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(comp_idx)%PTR%DOMAIN%&
             &TOPOLOGY%ELEMENTS
          DOMAIN_NODES=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES      
                    
          !write out the components' values of this node in this domain
          !DO scal_idx=1, NUM_OF_SCALING_FACTOR_SETS   
          IF(LIST_COMP_SCALE(comp_idx)==scal_idx) THEN
             scal_idx=scal_idx+1
             !global_number=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(elem_num)

            !finding the local numbering through the global to local mapping   
             !DOMAIN_MAPPING_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(LIST_COMP(scal_idx))%PTR%DOMAIN%MAPPINGS%ELEMENTS       
             !!get the domain index for this variable component according to my own computional node number
             !DO domain_idx=1,DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
             !   domain_no=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
             !   IF(domain_no==my_computational_node_number) THEN
             !      MY_DOMAIN_INDEX=domain_idx
             !      EXIT !out of loop--domain_idx
             !   ENDIF
             !ENDDO !domain_idX
             !local_number=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
             !!use local domain information find the out the maximum number of derivatives
             !DOMAIN_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_num)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%ELEMENTS
          
             ns=0
             IF(ALLOCATED(SCALING_BUFFER)) DEALLOCATE(SCALING_BUFFER)            
             ALLOCATE(SCALING_BUFFER(SUM(BASIS%NUMBER_OF_DERIVATIVES(1:BASIS%NUMBER_OF_NODES))),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate scale buffer in IO",ERR,ERROR,*999)
          
             BASIS=>DOMAIN_ELEMENTS%ELEMENTS(local_number)%BASIS   
                           
             !CALL DISTRIBUTED_VECTOR_DATA_GET(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(comp_idx)%PTR%FIELD%SCALINGS%SCALINGS(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS(comp_idx)%PTR% &
             !&SCALING_INDEX)%SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)                         
             IF(BASIS%DEGENERATE==.FALSE.) THEN
                 DO nn=1,BASIS%NUMBER_OF_NODES
                   np=DOMAIN_ELEMENTS%ELEMENTS(local_number)%ELEMENT_NODES(nn)
                   DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                      nk=DOMAIN_ELEMENTS%ELEMENTS(local_number)%ELEMENT_DERIVATIVES(mk,nn)
                      ny2=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk)
                      ns=ns+1
                      SCALING_BUFFER(ns)=1!SCALE_FACTORS(ny2)                      
                   ENDDO !mk
                ENDDO !nn
              ELSE
                CALL FLAG_ERROR("exporting degenerated nodes has not been implemented",ERR,ERROR,*999) 
             ENDIF  
             CALL FIELD_IO_FORTRAN_FILE_WRITE_DP(FILE_ID, SCALING_BUFFER, ns, ERR,ERROR,*999)
             !IF(ASSOCIATED(SCALING_BUFFER)) DEALLOCATE(SCALING_BUFFER)              
          ENDIF !(LIST_COMP_SCALE(comp_idx)==scal_idx)LIST_COMP(:)
       ENDDO !comp_idx=1         
    ENDDO!elem_idx    
    
    !close a file    
    CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)

    !release the temporary memory
    IF(ALLOCATED(SCALING_BUFFER)) DEALLOCATE(SCALING_BUFFER)
    IF(ALLOCATED(SCALE_FACTORS)) DEALLOCATE(SCALE_FACTORS)
    IF(ALLOCATED(NODAL_NUMBER)) DEALLOCATE(NODAL_NUMBER)
    IF(ALLOCATED(LIST_COMP_SCALE)) DEALLOCATE(LIST_COMP_SCALE)
        
    CALL EXITS("FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE")
    RETURN
999 CALL ERRORS("FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE",ERR,ERROR)
    CALL EXITS("FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE")
    RETURN 1  
  END SUBROUTINE FIELD_IO_EXPORT_ELEMENTS_INTO_LOCAL_FILE

  !
  !================================================================================================================================
  !  

  !>Sort the Elemental_info_set according to the type of field variable components 
  SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_SORT(LOCAL_PROCESS_ELEMENTAL_INFO_SET, my_computational_node_number, ERR,ERROR,*)      
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<elemental information in this process
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), ALLOCATABLE :: tmp_components(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_ELEMENTS !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS1, DOMAIN_ELEMENTS2! domain nodes
    INTEGER(INTG) :: domain_idx,domain_no, MY_DOMAIN_INDEX !temporary variable
    INTEGER(INTG) :: global_number1, local_number1, global_number2, local_number2
    INTEGER(INTG) :: component_idx, nn1, nn2 ! nn, tmp2, tmp1!temporary variable
    LOGICAL :: SWITCH

    !from now on, global numbering are used
    CALL ENTERS("FIELD_IO_ELEMENTAL_INFO_SET_SORT",ERR,ERROR,*999)    
  
    IF(.NOT.ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("list of global numbering in the input data is invalid",ERR,ERROR,*999)     
    ENDIF
    IF(.NOT.ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)) THEN
       CALL FLAG_ERROR("nodal information set in the input data is invalid",ERR,ERROR,*999)     
    ENDIF
 
    
    !!get my own computianal node number--be careful the rank of process in the MPI pool 
    !!is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
    !!will be a secured way to get the number         
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
    !IF(ERR/=0) GOTO 999        

    !group nodal information set according to its components, i.e. put all the nodes with the same components together
    !and change the global number in the LIST_OF_GLOBAL_NUMBER 
    nn1=1
    DO WHILE(nn1<LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS)
       !global number of this node
       global_number1=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn1)
       DO nn2=nn1+1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
          global_number2=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn2) 
          IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS==&
           &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS) THEN  
             SWITCH=.TRUE.
             !we will check the component (type of component, partial derivative).
             DO component_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS       
                !not safe, but it is fast            
            !=============================================================================================!
            !           checking according to local memory adddress                                       !
            !=============================================================================================!
            !are they in the same memory address?
            IF(.NOT.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR, &   
              &TARGET=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR))  THEN
               SWITCH=.FALSE. !out of loop-component_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS    
               EXIT
                ENDIF !ASSCOCIATED
            
                !! better use this one because it is safe method, but slow            
            !!=============================================================================================!
            !!           checking according to the types defined in the openCMISS                          !
            !!=============================================================================================!
            !!are they in the same field?
            !IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%FIELD%GLOBAL_NUMBER/= &
            !&LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%FIELD%GLOBAL_NUMBER) THEN
            !   SWITCH=.FALSE.
            !   EXIT
            !ELSE  !GLOBAL_NUBMER 
            !   !are they the same variable?
               !   IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%FIELD_VARIABLE%VARIABLE_NUMBER/= &
               !   & LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%FIELD_VARIABLE%VARIABLE_NUMBER) THEN
              !       SWITCH=.FALSE.
            !       EXIT
            !    ELSE !VARIABLE_NUBMER  
               !       !are they the same component?
              !      IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%COMPONENT_NUMBER/=&   
               !        &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%COMPONENT_NUMBER) THEN
              !          SWITCH=.FALSE.
             !          EXIT                   
            !       ENDIF !COMPONENT_NUMBER
            !   ENDIF ! VARIABLE_NUBMER
            !ENDIF !GLOBAL_NUBMER            
            ENDDO !component_idx
             
            !check whether correspoding two components have the same partial derivatives 
            IF(SWITCH==.TRUE.) THEN
               DO component_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS                    
               
                  !finding the local numbering for the NODAL_INFO_SET(nn1)   
                  DOMAIN_MAPPING_ELEMENTS=>&
                   &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%DOMAIN%MAPPINGS%ELEMENTS       
                   !get the domain index for this variable component according to my own computional node number
                  DO domain_idx=1,DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number1)%NUMBER_OF_DOMAINS
                     domain_no=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number1)%DOMAIN_NUMBER(domain_idx)
                     IF(domain_no==my_computational_node_number) THEN
                        MY_DOMAIN_INDEX=domain_idx
                        EXIT !out of loop--domain_idx
                     ENDIF
                  ENDDO !domain_idX
                  !local number of nn1'th node in the damain assoicated with component(component_idx)
                  local_number1=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number1)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
                  DOMAIN_ELEMENTS1=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%&
                      &DOMAIN%TOPOLOGY%ELEMENTS
               
                  !finding the local numbering for the NODAL_INFO_SET(nn2)   
                  DOMAIN_MAPPING_ELEMENTS=>&
                  &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%DOMAIN%MAPPINGS%ELEMENTS       
                  !get the domain index for this variable component according to my own computional node number
                  DO domain_idx=1,DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number2)%NUMBER_OF_DOMAINS
                     domain_no=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number2)%DOMAIN_NUMBER(domain_idx)
                     IF(domain_no==my_computational_node_number) THEN
                        MY_DOMAIN_INDEX=domain_idx
                        EXIT !out of loop--domain_idx
                     ENDIF
                  ENDDO !domain_idX
                  !local number of nn2'th node in the damain assoicated with component(component_idx)
                  local_number2=DOMAIN_MAPPING_ELEMENTS%GLOBAL_TO_LOCAL_MAP(global_number2)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
                  DOMAIN_ELEMENTS2=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%&
                      &DOMAIN%TOPOLOGY%ELEMENTS                   
                   
                  !checking whether they have the same basis
                  IF(DOMAIN_ELEMENTS1%ELEMENTS(local_number1)%BASIS%GLOBAL_NUMBER/=&
                     &DOMAIN_ELEMENTS2%ELEMENTS(local_number2)%BASIS%GLOBAL_NUMBER) THEN                    
                     SWITCH=.FALSE.
                     EXIT     
                  ENDIF   !DOMAIN_ELEMENTS1
               ENDDO !component_idx
            ENDIF !SWITCH==.TRUE.
         ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS==LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn+1)%NUMBER_OF_COMPONENTS
          
          !find two elements which have the same output, and then they should put together
          IF(SWITCH==.TRUE.) THEN
             !exchange the pointer(to the list the components)             
             IF(ALLOCATED(tmp_components)) DEALLOCATE(tmp_components)
             ALLOCATE(tmp_components(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)           
             tmp_components(:)=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(:)
             
             DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS)
             ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(&
                &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate buffer in elemental routine in FIELD IO",ERR,ERROR,*999)                         
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(:)=&
             &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%COMPONENTS(:)
             
             DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%COMPONENTS)
             ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%COMPONENTS(&
                &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate in elemental routine in FIELD IO",ERR,ERROR,*999)            
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%COMPONENTS(:)=tmp_components(:)
                          
             !setting the header information
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%SAME_HEADER=.FALSE.
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%SAME_HEADER=.TRUE.
             
             !exchange the number of components
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS=&
             &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS=&
             &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS

             !exchange the global number
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn2)=LOCAL_PROCESS_ELEMENTAL_INFO_SET%&
                 &LIST_OF_GLOBAL_NUMBER(nn1+1)
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn1+1)=global_number2
             
             !increase nn1 to skip the nodes which have the same output
             nn1=nn1+1
          ENDIF !(SWITCH=.TRUE.)                    
       ENDDO !nn2
       !increase the nn1 to check next node
       nn1=nn1+1        
    ENDDO !nn1<LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES   

    !order the variable components and group them: X1(1),X1(2),X1(3),X2(2),X2(3),X3(2)....
    !DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
    !   print "(A, I)", "nn=", nn
    !   !temporarily use nk, nu here to save memory
    !   IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=1) THEN
    !     component_idx=1      
    !     DO WHILE(component_idx<LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS) 
    !        !checking the same variable's components
    !       print "(A, I)", "component_idx=", component_idx
    !       print "(A, I)", "LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS", LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS
    !       DO WHILE(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(component_idx)%PTR%FIELD_VARIABLE, &
    !       & TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(component_idx+1)%PTR%FIELD_VARIABLE))
    !          component_idx=component_idx+1      
    !          IF(component_idx>=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS) THEN
    !             EXIT
    !          ENDIF                     
    !        ENDDO
    !       
    !        !It may have more than 3 component in the future?!! I do not know,too
    !        !so there the components are sorted according their numbering of component
    !        !nk and nu are used here temporarily
    !        DO tmp1=1,component_idx
    !           print "(A, I)", "tmp1=", tmp1
    !           SWITCH=.FALSE.
    !           DO tmp2=1,(component_idx-tmp1)
    !              IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2)%PTR%COMPONENT_NUMBER>&
    !              &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2+1)%PTR%COMPONENT_NUMBER) THEN                 
    !                 tmp_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2+1)%PTR
    !                 LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2+1)%PTR=>&
    !                 LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2)%PTR
    !                
    !                 LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2)%PTR=>tmp_ptr
    !                 SWITCH=.TRUE.
    !              ENDIF
    !           ENDDO
    !           IF(SWITCH==.TRUE.) THEN
    !              EXIT
    !           ENDIF  
    !        ENDDO
    !        NULLIFY(tmp_ptr)
    !        component_idx=component_idx+1
    !     ENDDO ! WHILE(component_idx<LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS)  
    !   ENDIF ! LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=1  
    !ENDDO !nn                 
    
    IF(ALLOCATED(tmp_components)) DEALLOCATE(tmp_components)
    
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_SORT")
    RETURN
999 CALL ERRORS("FIELD_IO_ELEMENTAL_INFO_SET_SORT",ERR,ERROR)
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_SORT")
    RETURN 1  
  END SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_SORT
  
  !
  !================================================================================================================================
  !  

  !>Collect the elemental information from each MPI process     
  SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), INTENT(INOUT):: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT):: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(DOMAIN_MAPPING_TYPE), POINTER:: DOMAIN_ELEMENTS_MAPPING !nodes in local mapping--it is different as exnode
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), ALLOCATABLE:: tmp_comp(:) !component of field variable
    TYPE(FIELD_VARIABLE_TYPE), POINTER:: FIELD_VARIABLE !field variable
    INTEGER(INTG), ALLOCATABLE:: NEW_LIST(:)
    INTEGER(INTG) :: field_idx, var_idx, component_idx, np, nn!temporary variable
    LOGICAL :: SWITCH

    CALL ENTERS("FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS",ERR,ERROR,*999)    
   
    !initialize the pointer
    !NULLIFY(NEW_LIST)
   
    !attache local process to local nodal information set. In current opencmiss system,
    !each local process owns it local nodal information, so all we need to do is to fill the nodal 
    !information set with nodal information of local process 
    IF((LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS==0).AND.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS) &
    & .AND.(.NOT.ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET))) THEN  
      DO field_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS%NUMBER_OF_FIELDS
         FIELD=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS%FIELDS(field_idx)%PTR
         DO var_idx=1, FIELD%NUMBER_OF_VARIABLES
          IF(ALLOCATED(FIELD%VARIABLES)) THEN           
               FIELD_VARIABLE=>FIELD%VARIABLES(var_idx)
            ELSE   
               EXIT
            ENDIF        
            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
               IF(ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%MAPPINGS%ELEMENTS
                  DO np=1,DOMAIN_ELEMENTS_MAPPING%NUMBER_OF_LOCAL
                     SWITCH=.FALSE.
                     DO nn=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
                        IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_ELEMENTS_MAPPING% &
                            &LOCAL_TO_GLOBAL_MAP(np)) THEN 
                           SWITCH=.TRUE.
                           EXIT
                        ENDIF                    
                     ENDDO
                     !have one more global node
                     !i hate the codes here, but i have to save the memory 
                     IF(SWITCH==.FALSE.) THEN
                        IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS/=0) THEN
                           !crease a new memory space
                           ALLOCATE(NEW_LIST(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS),STAT=ERR)
                           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                             !add one more node
                           NEW_LIST(1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS)=LOCAL_PROCESS_ELEMENTAL_INFO_SET%&
                               &LIST_OF_GLOBAL_NUMBER(1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS)              
                           !NEW_LIST(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS+1)=&
                           !&DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(np)
                           !release the old memory space
                           DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
                           ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_PROCESS_ELEMENTAL_INFO_SET%&
                           &NUMBER_OF_ELEMENTS+1),STAT=ERR)
                           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS)&
                           &=NEW_LIST(:)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS+1)&
                           &=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(np)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS+1
                           DEALLOCATE(NEW_LIST)                           
                        ELSE
                           ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1),STAT=ERR)
                           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1)=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(np)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=1
                        ENDIF   !LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES/=0)                 
                     ENDIF !SWITCH                                 
                 ENDDO !np
               ENDIF!ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES                             
            ENDDO !component_idx
         ENDDO !var_idx         
      ENDDO !field_idx
     
      !allocate the nodal information set and initialize them
      ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodal information set",ERR,ERROR,*999)    
      !ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES),STAT=ERR)
      !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodal information set",ERR,ERROR,*999)   
      !LOCAL_PROCESS_NODAL_INFO_SET%MAXIMUM_NUMBER_OF_DERIVATIVES=0 
      DO nn=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
         LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%SAME_HEADER=.FALSE.
         !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%LEN_OF_NODAL_INFO=0
         LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS=0
         IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS)) &
            &DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS)                                  
      ENDDO
   
      !collect nodal information from local process
      DO field_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS%NUMBER_OF_FIELDS
         FIELD=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS%FIELDS(field_idx)%PTR
         DO var_idx=1, FIELD%NUMBER_OF_VARIABLES
              IF(ALLOCATED(FIELD%VARIABLES)) THEN 
                FIELD_VARIABLE=>FIELD%VARIABLES(var_idx)
            ELSE   
               EXIT
            ENDIF        
            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
               IF(ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%MAPPINGS%ELEMENTS
                
                 DO np=1,DOMAIN_ELEMENTS_MAPPING%NUMBER_OF_LOCAL
                    DO nn=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
                       IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_ELEMENTS_MAPPING%&
                       &LOCAL_TO_GLOBAL_MAP(np)) THEN 
                          EXIT
                       ENDIF                    
                    ENDDO
                    !allocate variable component memory
                    IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=0) THEN
                        ALLOCATE(tmp_comp(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                        
                        tmp_comp(:)=LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(:)
                        DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS)
                        ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%&
                           &COMPONENTS(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                        LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(1:&
                        &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS)=tmp_comp(:)
                        LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(&
                        &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1)%PTR=>&
                        FIELD_VARIABLE%COMPONENTS(component_idx)
                        DEALLOCATE(tmp_comp)
                        
                        !tmp_comp=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS
                        !NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS)
                        !LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(&
                        !   &1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS)=tmp_comp(:)
                 
                        !LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(LOCAL_PROCESS_ELEMENTAL_INFO_SET%&
                        !   &ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1)%PTR=>                                       
                     ELSE
                        ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(1),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                        LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(1)%PTR=>&
                            &FIELD_VARIABLE%COMPONENTS(component_idx)
                     ENDIF
                     !increase number of component
                     LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS=&
                     &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1                    
                  ENDDO !np                       
               ENDIF   !(ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES))          
              
               !we do not need to check the nodes
               !DO np=1,DOMAIN_NODES%NUMBER_OF_NODES                 
               !   DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
               !      IF(LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER) THEN                    
               !         EXIT
               !      ENDIF                                        
               !   ENDDO                 
               !ENDDO !np
            ENDDO !component_idx
        ENDDO !var_idx         
      ENDDO !field_idx     
      !NULLIFY(tmp_comp)
     
      !LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>LIST_OF_GLOBAL_NUMBER
      !NULLIFY(LIST_OF_GLOBAL_NUMBER)      
      !release the temporary meomery           
      IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)             
    ELSE
      CALL FLAG_ERROR("nodal information set is not initialized properly, and call start method first",ERR,ERROR,*999)
    ENDIF
    
    IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    IF(ALLOCATED(tmp_comp)) DEALLOCATE(tmp_comp)
    
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN
999 CALL ERRORS("FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS",ERR,ERROR)
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN 1  
  END SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS      

  !
  !================================================================================================================================
  !  

  !>Finalized the elemental information set     
  SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nele, ncomp  !temporary variable

    CALL ENTERS("FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE",ERR,ERROR,*999)    
   
    IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)) THEN
       NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)
    ENDIF
    DO nele=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
       IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)) THEN
          DO ncomp=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%NUMBER_OF_COMPONENTS
             IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS(ncomp)%PTR)) THEN
                NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS(ncomp)%PTR)
             ENDIF
          ENDDO
          LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%NUMBER_OF_COMPONENTS=0
          LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%SAME_HEADER=.FALSE.
          IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)) THEN
             DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)
          ENDIF   
       ENDIF
    ENDDO
   
    LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=0
    IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    ENDIF
            
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE")
    RETURN
999 CALL ERRORS("FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE",ERR,ERROR)
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE")
    RETURN 1  
  END SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_FINALIZE

  !
  !================================================================================================================================
  !  

  !>Initialize the elemental information set         
  SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FIELDS, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<elemental information in this process
    TYPE(FIELDS_TYPE), POINTER ::FIELDS !<the field object
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: nfid, nele, ncomp  !temporary variable

    CALL ENTERS("FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE",ERR,ERROR,*999)    

    !validate the input data
    !checking whether the list of fields in the same region
    IF(ASSOCIATED(FIELDS%REGION)) THEN
       DO nfid =2, FIELDS%NUMBER_OF_FIELDS
          IF(FIELDS%FIELDS(nfid-1)%PTR%REGION%USER_NUMBER/=FIELDS%FIELDS(nfid)%PTR%REGION%USER_NUMBER) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(nfid-1,"*",ERR,ERROR))//" and "//&
                 &TRIM(NUMBER_TO_VSTRING(nfid,"*",ERR,ERROR))//" fields are not in the same region"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)  
             GOTO 999       
          ENDIF       
       ENDDO 
     
       !checking whether the field's handle 
       DO nfid=1, FIELDS%NUMBER_OF_FIELDS
          IF(.NOT.ASSOCIATED(FIELDS%FIELDS(nfid)%PTR)) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(nfid,"*",ERR,ERROR))//" field handle in fields list is invalid"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)  
             GOTO 999       
          ENDIF       
       ENDDO
       !checking whether the list of fields are using the same decomposition     
       !IF(.NOT.ASSOCIATED(DECOMPOSITION))
       !  CALL FLAG_ERROR("decomposition method is not vakid",ERR,ERROR,*999)
       !ENDIF
       !DO nfid =1, FIELDS%NUMBER_OF_FIELDS
       !  IF(FIELDS%FIELDS(nfid)%PTR%DECOMPOSITION/=DECOMPOSITION)
       !    LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(nfid,"*",ERR,ERROR)) //" field "&
       !    & //" uses different decomposition method with the specified decomposition method,"//&
       !     & "which is not supported currently, ask Heye for more details"       
       !    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
       !  ENDIF       
       !ENDDO                   
    ELSE
      CALL FLAG_ERROR("list of Field is not associated with any region",ERR,ERROR,*999)
      GOTO 999
    ENDIF !ASSOCIATED(FIELDS%REGION
   
    !release the pointers if they are associated
    !Allocatable is the refer to dynamic size of array
    !IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET)) THEN
       IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)) THEN
          DO nele=1, LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
             DO ncomp=1, LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%NUMBER_OF_COMPONENTS
                NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS(ncomp)%PTR)
             ENDDO   
             DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)
          ENDDO
          DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)
       ENDIF
    !ELSE
    !   ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET,STAT=ERR)
    !   IF(ERR/=0) CALL FLAG_ERROR("could not allocate nodal information set for IO",ERR,ERROR,*999)             
    !ENDIF
   
    !set to number of nodes to zero
    LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=0 
    IF(ALLOCATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) &
    & DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    
    !associated nodel info set with the list of fields
    IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)) THEN
       NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)
    ENDIF
    LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS=>FIELDS
 
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE")
    RETURN 1  
  END SUBROUTINE FIELD_IO_ELEMENTAL_INFO_SET_INITIALISE

  !!
  !!================================================================================================================================
  !!
  !
  !!>Import nodal information \see{FIELD_IO::FIELD_IO_NODES_IMPORT}.                 
  !SUBROUTINE FIELD_IO_NODES_IMPORT(FIELDS, FILE_NAME, METHOD, ERR,ERROR,*)
  !  !Argument variables       
  !  TYPE(FIELDS_TYPE), POINTER :: FIELDS !<the field object
  !  TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<file name
  !  TYPE(VARYING_STRING), INTENT(IN):: METHOD
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(FIELD_IO_NODAL_INFO_SET) :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
  !  INTEGER(INTG):: my_computational_node_number !<local process number    
  !  INTEGER(INTG):: computational_node_numbers   !<total process number      
  !
  !  CALL ENTERS("FIELD_IO_NODES_IMPORT", ERR,ERROR,*999)    
  !
  !  !Get the number of computational nodes
  !  computational_node_numbers=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  !  IF(ERR/=0) GOTO 999
  !  !Get my computational node number
  !  my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  !  IF(ERR/=0) GOTO 999
  !  IF(METHOD=="FORTRAN") THEN
  !     CALL FIELD_IO_NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_NODAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
  !     CALL FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)
  !     CALL FIELD_IO_NODAL_INFO_SET_SORT(LOCAL_PROCESS_NODAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)    
  !     CALL FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE(LOCAL_PROCESS_NODAL_INFO_SET, FILE_NAME, my_computational_node_number, &
  !          &computational_node_numbers, ERR, ERROR, *999)
  !     CALL FIELD_IO_NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)
  !  ELSE IF(METHOD=="MPIIO") THEN
  !     CALL FLAG_ERROR("what are u thinking, of course not!",ERR,ERROR,*999)
  !     !CALL FIELD_IO_NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_NODAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
  !     !CALL FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)
  !     !CALL FIELD_IO_NODAL_INFO_SET_SORT(LOCAL_PROCESS_NODAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)
  !     !CALL FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE(LOCAL_PROCESS_NODAL_INFO_SET, FILE_NAME, my_computational_node_number, &
  !     !     &computational_node_numbers, ERR, ERROR, *999)
  !     !CALL FIELD_IO_NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)           
  !  ENDIF
  !  
  !  CALL EXITS("FIELD_IO_NODES_IMPORT")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_NODES_IMPORT",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_NODES_IMPORT")
  !  RETURN 1  
  !END SUBROUTINE FIELD_IO_NODES_IMPORT


  !
  !================================================================================================================================
  !

  !>Export nodal information \see{FIELD_IO::FIELD_IO_NODES_EXPORT}.                 
  SUBROUTINE FIELD_IO_NODES_EXPORT(FIELDS, FILE_NAME, METHOD, ERR,ERROR,*)
    !Argument variables       
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<the field object
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<file name
    TYPE(VARYING_STRING), INTENT(IN):: METHOD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_IO_NODAL_INFO_SET) :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG):: my_computational_node_number !<local process number    
    INTEGER(INTG):: computational_node_numbers   !<total process number      

    CALL ENTERS("FIELD_IO_NODES_EXPORT", ERR,ERROR,*999)    

    !Get the number of computational nodes
    computational_node_numbers=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !Get my computational node number
    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    IF(METHOD=="FORTRAN") THEN
       CALL FIELD_IO_NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_NODAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
       CALL FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)
       CALL FIELD_IO_NODAL_INFO_SET_SORT(LOCAL_PROCESS_NODAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)    
       CALL FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE(LOCAL_PROCESS_NODAL_INFO_SET, FILE_NAME, my_computational_node_number, &
            &computational_node_numbers, ERR, ERROR, *999)
       CALL FIELD_IO_NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)
    ELSE IF(METHOD=="MPIIO") THEN
       CALL FLAG_ERROR("MPI IO has not been implemented yet!",ERR,ERROR,*999)
       !CALL FIELD_IO_NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_NODAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
       !CALL FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)
       !CALL FIELD_IO_NODAL_INFO_SET_SORT(LOCAL_PROCESS_NODAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)
       !CALL FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE(LOCAL_PROCESS_NODAL_INFO_SET, FILE_NAME, my_computational_node_number, &
       !     &computational_node_numbers, ERR, ERROR, *999)
       !CALL FIELD_IO_NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*999)           
    ELSE 
       CALL FLAG_ERROR("Unknown method!",ERR,ERROR,*999)   
    ENDIF   
    
    CALL EXITS("FIELD_IO_NODES_EXPORT")
    RETURN
999 CALL ERRORS("FIELD_IO_NODES_EXPORT",ERR,ERROR)
    CALL EXITS("FIELD_IO_NODES_EXPORT")
    RETURN 1  
  END SUBROUTINE FIELD_IO_NODES_EXPORT

  !
  !================================================================================================================================
  !  

  !>Initialize nodal information set 
  SUBROUTINE FIELD_IO_NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_NODAL_INFO_SET, FIELDS, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    TYPE(FIELDS_TYPE), POINTER ::FIELDS !<the field object
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: num_field, nn, ncomp  !temporary variable

    CALL ENTERS("FIELD_IO_NODAL_INFO_SET_INITIALISE",ERR,ERROR,*999)    

    !validate the input data
    !checking whether the list of fields in the same region
    IF(ASSOCIATED(FIELDS%REGION)) THEN
       DO num_field =2, FIELDS%NUMBER_OF_FIELDS
          IF(FIELDS%FIELDS(num_field-1)%PTR%REGION%USER_NUMBER/=FIELDS%FIELDS(num_field)%PTR%REGION%USER_NUMBER) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(num_field-1,"*",ERR,ERROR))//" and "//&
             &TRIM(NUMBER_TO_VSTRING(num_field,"*",ERR,ERROR))//" fields are not in the same region"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)     
             GOTO 999    
          ENDIF       
       ENDDO 
       
       !checking whether the field's handle 
       DO num_field=1, FIELDS%NUMBER_OF_FIELDS
          IF(.NOT.ASSOCIATED(FIELDS%FIELDS(num_field)%PTR)) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(num_field,"*",ERR,ERROR))//" field handle in fields list is invalid"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
             GOTO 999
          ENDIF       
       ENDDO
       !checking whether the list of fields are using the same decomposition     
       !IF(.NOT.ASSOCIATED(DECOMPOSITION))
       !  CALL FLAG_ERROR("decomposition method is not vakid",ERR,ERROR,*999)
       !ENDIF
       !DO num_field =1, FIELDS%NUMBER_OF_FIELDS
       !  IF(FIELDS%FIELDS(num_field)%PTR%DECOMPOSITION/=DECOMPOSITION)
       !    LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(num_field,"*",ERR,ERROR)) //" field "&
       !    & //" uses different decomposition method with the specified decomposition method,"//&
       !     & "which is not supported currently, ask Heye for more details"       
       !    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
       !  ENDIF       
       !ENDDO                   
    ELSE
       CALL FLAG_ERROR("list of Field is not associated with any region",ERR,ERROR,*999)
       GOTO 999
    ENDIF !ASSOCIATED(FIELDS%REGION
   
    !release the pointers if they are associated
    !Allocatable is the refer to dynamic size of arrayNode:
    !IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET)) THEN
       IF(ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
          DO nn=1, LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
             DO ncomp=1, LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS
               NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(ncomp)%PTR)
             ENDDO   
             DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
          ENDDO
          DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)
       ENDIF
    !ELSE
    !   ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET,STAT=ERR)
    !   IF(ERR/=0) CALL FLAG_ERROR("could not allocate nodal information set for IO",ERR,ERROR,*999)             
    !ENDIF
   
    !set to number of nodes to zero
    LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=0 
    IF(ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    
    !associated nodel info set with the list of fields
    IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)) NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)
    LOCAL_PROCESS_NODAL_INFO_SET%FIELDS=>FIELDS
 
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_IO_NODAL_INFO_SET_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_INITIALISE")
    RETURN 1  
  END SUBROUTINE FIELD_IO_NODAL_INFO_SET_INITIALISE

  !
  !================================================================================================================================
  !  

  !>Collect nodal information from each MPI process      
  SUBROUTINE FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT):: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT):: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(DOMAIN_NODES_TYPE), POINTER:: DOMAIN_NODES !nodes in local domain
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), ALLOCATABLE:: tmp_comp(:) !component of field variable
    TYPE(FIELD_VARIABLE_TYPE), POINTER:: FIELD_VARIABLE !field variable
    INTEGER(INTG) :: field_idx, var_idx, component_idx, np, nn!temporary variable
    INTEGER(INTG), ALLOCATABLE:: NEW_LIST(:)
    LOGICAL :: SWITCH

    CALL ENTERS("FIELD_IO_NODAL_INFO_SET_CREATE_FINISH",ERR,ERROR,*999)    
   
    !initialize the pointer
    !NULLIFY(NEW_LIST)
   
    !attache local process to local nodal information set. In current opencmiss system,
    !each local process owns it local nodal information, so all we need to do is to fill the nodal 
    !information set with nodal information of local process 
    IF((LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0).AND.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS) &
    & .AND.(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET))) THEN  
       DO field_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%FIELDS%NUMBER_OF_FIELDS
          FIELD=>LOCAL_PROCESS_NODAL_INFO_SET%FIELDS%FIELDS(field_idx)%PTR
          DO var_idx=1, FIELD%NUMBER_OF_VARIABLES
             IF(ALLOCATED(FIELD%VARIABLES)) THEN           
                FIELD_VARIABLE=>FIELD%VARIABLES(var_idx)
             ELSE   
                EXIT
             ENDIF        
             DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES)) THEN
                   DOMAIN_NODES=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES
                   DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
                      SWITCH=.FALSE.
                      DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
                         IF(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER) THEN
                            SWITCH=.TRUE.
                            EXIT
                         ENDIF                    
                      ENDDO!nn=1
                      !have one more global node
                      !i hate the codes here, but i have to save the memory 
                      IF(SWITCH==.FALSE.) THEN
                         IF(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES/=0) THEN
                            !crease a new memory space
                            ALLOCATE(NEW_LIST(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                            !add one more node
                            NEW_LIST(:)=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(:)
                            DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
                            ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_PROCESS_NODAL_INFO_SET% &
                              &NUMBER_OF_NODES+1),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                            LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1:LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES)=&
                            &NEW_LIST(:)
                            LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1)=&
                            &DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER                            
                            DEALLOCATE(NEW_LIST)
                            
                            !NEW_LIST(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1)=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
                            !!release the old memory space
                            !DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
                            !LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>NEW_LIST
                            !NULLIFY(NEW_LIST)
                            !LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1
                         ELSE
                            ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                            LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1)=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
                         ENDIF   !LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES/=0)   
                         LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1                                       
                      ENDIF !SWITCH                                 
                   ENDDO !np
                ENDIF!ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES                             
             ENDDO !component_idx
          ENDDO !var_idx         
       ENDDO !field_idx
     
       !allocate the nodal information set and initialize them
       ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES),STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodal information set",ERR,ERROR,*999)    
       !ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES),STAT=ERR)
       !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodal information set",ERR,ERROR,*999)   
       !LOCAL_PROCESS_NODAL_INFO_SET%MAXIMUM_NUMBER_OF_DERIVATIVES=0 
       DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
          LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER=.FALSE.
          !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%LEN_OF_NODAL_INFO=0
          LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS=0
          IF(ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)) &
             &DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)                                  
       ENDDO
   
       !collect nodal information from local process
       DO field_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%FIELDS%NUMBER_OF_FIELDS
          FIELD=>LOCAL_PROCESS_NODAL_INFO_SET%FIELDS%FIELDS(field_idx)%PTR
          DO var_idx=1, FIELD%NUMBER_OF_VARIABLES
             IF(ALLOCATED(FIELD%VARIABLES)) THEN 
                FIELD_VARIABLE=>FIELD%VARIABLES(var_idx)
             ELSE   
                EXIT
             ENDIF        
             DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES)) THEN
                   DOMAIN_NODES=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES
               
                   DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
                      DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
                         IF(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER) THEN
                            EXIT
                         ENDIF                    
                      ENDDO
                      !allocate variable component memory
                      IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=0) THEN
                         !crease a new memory space
                         ALLOCATE(tmp_comp(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS),STAT=ERR)
                         IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                         !add one more node
                         tmp_comp(:)=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(:)
                         DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
                         ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(&
                            &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1),STAT=ERR)
                         IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                         LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(1:&
                            &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS)=tmp_comp(:)
                         LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)&
                            &%NUMBER_OF_COMPONENTS+1)%PTR=>FIELD_VARIABLE%COMPONENTS(component_idx)
                         DEALLOCATE(tmp_comp)
                            
                         !tmp_comp=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS
                         !NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
                      
                         !ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(&
                         ! &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1),STAT=ERR)
                         !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                
                         !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(1:LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)&
                         !&%NUMBER_OF_COMPONENTS)=tmp_comp(:)

                         !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)&
                         !&%NUMBER_OF_COMPONENTS+1)%PTR=>FIELD_VARIABLE%COMPONENTS(component_idx)                                       
                      ELSE
                         ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(1),STAT=ERR)
                         IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                         LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(1)%PTR=>FIELD_VARIABLE%COMPONENTS(component_idx)
                      ENDIF
                      !increase number of component
                      LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS=&
                      &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1                    
                   ENDDO !np                       
                ENDIF   !(ASSOCIATED(FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES))          
             
                !we do not need to check the nodes
                !DO np=1,DOMAIN_NODES%NUMBER_OF_NODES                 
                !   DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
                !      IF(LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER) THEN                    
                !         EXIT
                !      ENDIF                                        
                !   ENDDO                 
                !ENDDO !np
             ENDDO !component_idx
          ENDDO !var_idx         
       ENDDO !field_idx     
       !NULLIFY(tmp_comp)
     
       !LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>LIST_OF_GLOBAL_NUMBER
       !NULLIFY(LIST_OF_GLOBAL_NUMBER)      
       !release the temporary meomery           
       IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
       IF(ALLOCATED(tmp_comp)) DEALLOCATE(NEW_LIST)
    ELSE
     CALL FLAG_ERROR("nodal information set is not initialized properly, and call start method first",ERR,ERROR,*999)
    ENDIF !(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0).AND.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)

    IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    IF(ALLOCATED(tmp_comp)) DEALLOCATE(NEW_LIST)
 
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN
999 CALL ERRORS("FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS",ERR,ERROR)
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN 1  
  END SUBROUTINE FIELD_IO_NODAL_INFO_SET_ATTACH_LOCAL_PROCESS      

  !
  !================================================================================================================================
  !  

  !>Sort nodal information according to the type of field variable component           
  SUBROUTINE FIELD_IO_NODAL_INFO_SET_SORT(LOCAL_PROCESS_NODAL_INFO_SET, my_computational_node_number, ERR,ERROR,*)      
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), ALLOCATABLE :: tmp_components(:)        
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES1, DOMAIN_NODES2! domain nodes
    INTEGER(INTG) :: domain_idx,domain_no, MY_DOMAIN_INDEX !temporary variable
    INTEGER(INTG) :: global_number1, local_number1, global_number2, local_number2
    INTEGER(INTG) :: component_idx, tmp1, nn1, nn2 ! nn, tmp2!temporary variable
    INTEGER(INTG), ALLOCATABLE:: array1(:), array2(:)  
    LOGICAL :: SWITCH

    !from now on, global numbering are used
    CALL ENTERS("FIELD_IO_NODAL_INFO_SET_SORT",ERR,ERROR,*999)    
  
    IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("list of global numbering in the input data is invalid",ERR,ERROR,*999)     
    ENDIF
    IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
       CALL FLAG_ERROR("nodal information set in the input data is invalid",ERR,ERROR,*999)     
    ENDIF
 
    
    !!get my own computianal node number--be careful the rank of process in the MPI pool 
    !!is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
    !!will be a secured way to get the number         
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
    !IF(ERR/=0) GOTO 999        

    !group nodal information set according to its components, i.e. put all the nodes with the same components together
    !and change the global number in the LIST_OF_GLOBAL_NUMBER 
    nn1=1
    DO WHILE(nn1<LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES)
       !global number of this node
       global_number1=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn1)
       DO nn2=nn1+1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
          global_number2=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn2) 
          IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS==&
           &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS) THEN  
             SWITCH=.TRUE.
             !we will check the component (type of component, partial derivative).
             DO component_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS       
                !!not safe, but it is fast            
            !!=============================================================================================!
            !!           checking according to local memory adddress                                       !
            !!=============================================================================================!
            !!are they in the same memory address?
            !IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR/=&   
            !  &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR) 
            !THEN
            !   SWITCH=.FALSE. !out of loop-component_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS    
            !   EXIT
                !ENDIF                      NUMBER_OF_NODES
            
                ! better use this one because it is safe method, but slow            
            !=============================================================================================!
            !           checking according to the types defined in the openCMISS                          !
            !=============================================================================================!
            !are they in the same field?
            IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%FIELD%GLOBAL_NUMBER/= &
            &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%FIELD%GLOBAL_NUMBER) THEN
               SWITCH=.FALSE.
               EXIT
            ELSE  !GLOBAL_NUBMER 
               !are they the same variable?
                  IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%FIELD_VARIABLE%VARIABLE_NUMBER/= &
                  & LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%FIELD_VARIABLE%VARIABLE_NUMBER) THEN
                     SWITCH=.FALSE.
                   EXIT
                ELSE !VARIABLE_NUBMER  
                     !are they the same component?
                    IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%COMPONENT_NUMBER/=&   
                       &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%COMPONENT_NUMBER) THEN
                        SWITCH=.FALSE.
                       EXIT                   
                   ENDIF !COMPONENT_NUMBER
               ENDIF ! VARIABLE_NUBMER
            ENDIF !GLOBAL_NUBMER            
             ENDDO !component_idx
             
             !check whether correspoding two components have the same partial derivatives 
             IF(SWITCH==.TRUE.) THEN
                DO component_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS                    
               
                  !finding the local numbering for the NODAL_INFO_SET(nn1)   
                   DOMAIN_MAPPING_NODES=>&
                   &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%DOMAIN%MAPPINGS%NODES       
                   !get the domain index for this variable component according to my own computional node number
                   DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number1)%NUMBER_OF_DOMAINS
                      domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number1)%DOMAIN_NUMBER(domain_idx)
                      IF(domain_no==my_computational_node_number) THEN
                         MY_DOMAIN_INDEX=domain_idx
                         EXIT !out of loop--domain_idx
                      ENDIF
                   ENDDO !domain_idX
                   !local number of nn1'th node in the damain assoicated with component(component_idx)
                   local_number1=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number1)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
                   DOMAIN_NODES1=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%DOMAIN%TOPOLOGY%NODES
               
               !finding the local numbering for the NODAL_INFO_SET(nn2)   
                   DOMAIN_MAPPING_NODES=>&
                   &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%DOMAIN%MAPPINGS%NODES       
                   !get the domain index for this variable component according to my own computional node number
                   DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number2)%NUMBER_OF_DOMAINS
                      domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number2)%DOMAIN_NUMBER(domain_idx)
                      IF(domain_no==my_computational_node_number) THEN
                         MY_DOMAIN_INDEX=domain_idx
                         EXIT !out of loop--domain_idx
                      ENDIF
                   ENDDO !domain_idX
                   !local number of nn2'th node in the damain assoicated with component(component_idx)
                   local_number2=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number2)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
                   DOMAIN_NODES2=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%DOMAIN%TOPOLOGY%NODES                   
                   
                   !checking whether they have the same number of partiabl derivative
                   IF(DOMAIN_NODES1%NODES(local_number1)%NUMBER_OF_DERIVATIVES&
                      &==DOMAIN_NODES2%NODES(local_number2)%NUMBER_OF_DERIVATIVES) THEN                    
                 ALLOCATE(array1(DOMAIN_NODES1%NODES(local_number1)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                 IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty buffer in IO sorting",ERR,ERROR,*999)                                    
                 ALLOCATE(array2(DOMAIN_NODES1%NODES(local_number2)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                 IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty buffer in IO sorting",ERR,ERROR,*999)                                    
                      array1(1:DOMAIN_NODES1%NODES(local_number1)%NUMBER_OF_DERIVATIVES)=0
                      array2(1:DOMAIN_NODES1%NODES(local_number2)%NUMBER_OF_DERIVATIVES)=0                      
                      array1(1:DOMAIN_NODES1%NODES(local_number1)%NUMBER_OF_DERIVATIVES)=&
                      &DOMAIN_NODES1%NODES(local_number1)%PARTIAL_DERIVATIVE_INDEX(:)
                      array2(1:DOMAIN_NODES1%NODES(local_number2)%NUMBER_OF_DERIVATIVES)=&
                      &DOMAIN_NODES1%NODES(local_number2)%PARTIAL_DERIVATIVE_INDEX(:)
                    CALL LIST_SORT(array1,ERR,ERROR,*999)
                    CALL LIST_SORT(array2,ERR,ERROR,*999)  
                    tmp1=SUM(array1-array2)
                    DEALLOCATE(array1)
                    DEALLOCATE(array2)                    
                    IF(tmp1/=0) THEN
                         SWITCH=.FALSE.
                         EXIT !out of loop-component_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS    
                    ENDIF                  
                   ELSE
                      SWITCH=.FALSE.
                      EXIT !out of loop-component_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS    
                   ENDIF                                     
                ENDDO !component_idx
             ENDIF !SWITCH==.TRUE.
          ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS==LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn+1)%NUMBER_OF_COMPONENTS
          
          !find two nodes which have the same output, and then they should put together
          IF(SWITCH==.TRUE.) THEN
             !exchange the pointer(to the list the components)
             IF(ALLOCATED(tmp_components)) DEALLOCATE(tmp_components)
             ALLOCATE(tmp_components(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
             tmp_components(:)=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(:)
             
             DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS)
             ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(&
                &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS(:)=&
                &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%COMPONENTS(:)

             DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%COMPONENTS)
             ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%COMPONENTS(&
                &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)             
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%COMPONENTS(:)=tmp_components(:)             
              
             !setting the header information
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%SAME_HEADER=.FALSE.
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%SAME_HEADER=.TRUE.
             
             !exchange the number of components
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS=&
             &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS=&
             &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS

             !exchange the global number
             LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn2)=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn1+1)
             LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn1+1)=global_number2
             
             !increase nn1 to skip the nodes which have the same output
             nn1=nn1+1
          ENDIF !(SWITCH=.TRUE.)   
                 
       ENDDO !nn2
       !increase the nn1 to check next node
       nn1=nn1+1        
    ENDDO !nn1<LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES   

    !order the variable components and group them: X1(1),X1(2),X1(3),X2(2),X2(3),X3(2)....
    !DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
    !   print "(A, I)", "nn=", nn
    !   !temporarily use nk, nu here to save memory
    !   IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=1) THEN
    !     component_idx=1      
    !     DO WHILE(component_idx<LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONNode:ENTS) 
    !        !checking the same variable's components
    !       print "(A, I)", "component_idx=", component_idx
    !       print "(A, I)", "LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS", LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS
    !       DO WHILE(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(component_idx)%PTR%FIELD_VARIABLE, &
    !       & TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(component_idx+1)%PTR%FIELD_VARIABLE))
    !          component_idx=component_idx+1      
    !          IF(component_idx>=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS) THEN
    !             EXIT
    !          ENDIF                     
    !        ENDDO
    !       
    !        !It may have more than 3 component in the future?!! I do not know,too
    !        !so there the components are sorted according their numbering of component
    !        !nk and nu are used here temporarily
    !        DO tmp1=1,component_idx
    !           print "(A, I)", "tmp1=", tmp1
    !           SWITCH=.FALSE.
    !           DO tmp2=1,(component_idx-tmp1)
    !              IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2)%PTR%COMPONENT_NUMBER>&
    !              &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2+1)%PTR%COMPONENT_NUMBER) THEN                 
    !                 tmp_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2+1)%PTR
    !                 LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2+1)%PTR=>&
    !                 LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2)%PTR
    !                
    !                 LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(tmp2)%PTR=>tmp_ptr
    !                 SWITCH=.TRUE.
    !              ENDIF
    !           ENDDO
    !           IF(SWITCH==.TRUE.) THEN
    !              EXIT
    !           ENDIF  
    !        ENDDO
    !        NULLIFY(tmp_ptr)
    !        component_idx=component_idx+1
    !     ENDDO ! WHILE(component_idx<LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS)  
    !   ENDIF ! LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=1  
    !ENDDO !nn                 
    
    IF(ALLOCATED(tmp_components)) DEALLOCATE(tmp_components)
    IF(ALLOCATED(array1)) DEALLOCATE(array1)
    IF(ALLOCATED(array2)) DEALLOCATE(array2)

    CALL EXITS("FIELD_IO_NODAL_INFO_SET_SORT")
    RETURN
999 CALL ERRORS("FIELD_IO_NODAL_INFO_SET_SORT",ERR,ERROR)
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_SORT")
    RETURN 1  
  END SUBROUTINE FIELD_IO_NODAL_INFO_SET_SORT

  !
  !================================================================================================================================
  !  

  !>Finalize nodal information set    
  SUBROUTINE FIELD_IO_NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn, nm  !temporary variable

    CALL ENTERS("FIELD_IO_NODAL_INFO_SET_FINALIZE",ERR,ERROR,*999)    
   
    IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)) THEN
       NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)
    ENDIF
    DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
       IF(ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)) THEN
          DO nm=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS
             IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(nm)%PTR)) THEN
                NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(nm)%PTR)
             ENDIF
          ENDDO
          LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS=0
          LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER=.FALSE.
          IF(ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)) THEN
             DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
          ENDIF   
       ENDIF
    ENDDO
   
    LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=0
    IF(ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    ENDIF
            
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_FINALIZE")
    RETURN
999 CALL ERRORS("FIELD_IO_NODAL_INFO_SET_FINALIZE",ERR,ERROR)
    CALL EXITS("FIELD_IO_NODAL_INFO_SET_FINALIZE")
    RETURN 1  
  END SUBROUTINE FIELD_IO_NODAL_INFO_SET_FINALIZE

  !
  !================================================================================================================================
  !  

  !>Get the derivative information               
  FUNCTION FIELD_IO_LABEL_DERIVATIVE_INFO_GET(GROUP_DERIVATIVES, NUMBER_DERIVATIVES, LABEL_TYPE, ERR, ERROR)
    !Argument variables   
    INTEGER(INTG), INTENT(IN) :: NUMBER_DERIVATIVES
    INTEGER(INTG), INTENT(IN) :: GROUP_DERIVATIVES(NUMBER_DERIVATIVES)
    INTEGER(INTG), INTENT(IN) :: LABEL_TYPE !<identitor for information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) ::FIELD_IO_LABEL_DERIVATIVE_INFO_GET
    INTEGER(INTG) :: dev_idx
   
    CALL ENTERS("FIELD_IO_LABEL_DERIVATIVE_INFO_GET",ERR,ERROR,*999)    
   
    IF(NUMBER_DERIVATIVES==0) THEN
       CALL FLAG_ERROR("number of derivatives in the input data is zero",ERR,ERROR,*999)           
    ENDIF
    IF(LABEL_TYPE/=FIELD_IO_DERIVATIVE_LABEL) THEN
       CALL FLAG_ERROR("label type in the input data is not derivative label",ERR,ERROR,*999)           
    ENDIF

    IF((NUMBER_DERIVATIVES==1).AND.GROUP_DERIVATIVES(1)==NO_PART_DERIV) THEN
       FIELD_IO_LABEL_DERIVATIVE_INFO_GET=""
    ELSE
       FIELD_IO_LABEL_DERIVATIVE_INFO_GET="("
       DO dev_idx=1,NUMBER_DERIVATIVES
          SELECT CASE(GROUP_DERIVATIVES(dev_idx))    
            CASE(NO_PART_DERIV)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET
            CASE(PART_DERIV_S1)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d/ds1"      
            CASE(PART_DERIV_S1_S1)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds1ds1"       
            CASE(PART_DERIV_S2)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d/ds2"
            CASE(PART_DERIV_S2_S2)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds2ds2"
            CASE(PART_DERIV_S1_S2)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d/ds3"
            CASE(PART_DERIV_S3)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds3ds3"
            CASE(PART_DERIV_S3_S3)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds3ds3"
            CASE(PART_DERIV_S1_S3)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds1ds3"
            CASE(PART_DERIV_S2_S3)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds2ds3"
            CASE(PART_DERIV_S1_S2_S3)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d3/ds1ds2ds3"
            CASE(PART_DERIV_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d/ds4"
            CASE(PART_DERIV_S4_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds4ds4"
            CASE(PART_DERIV_S1_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds1ds4"
            CASE(PART_DERIV_S2_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds2ds4"
            CASE(PART_DERIV_S3_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d2/ds3ds4"
            CASE(PART_DERIV_S1_S2_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d3/ds1ds2ds4"
            CASE(PART_DERIV_S1_S3_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d3/ds1ds3ds4"
            CASE(PART_DERIV_S2_S3_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d3/ds2ds3ds4"
            CASE(PART_DERIV_S1_S2_S3_S4)
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET=FIELD_IO_LABEL_DERIVATIVE_INFO_GET//", d4/ds1ds2ds3ds4"
            CASE DEFAULT
              FIELD_IO_LABEL_DERIVATIVE_INFO_GET="unknown field variable type, add more details later, #Components="!&
              !&//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
         END SELECT
       ENDDO! dev_idx
    ENDIF !NUMBER_DERIVATIVES==1.AND.GROUP_DERIVATIVES(1)==NO_PART_DERIV
             
    CALL EXITS("FIELD_IO_LABEL_DERIVATIVE_INFO_GET")
    RETURN
999 CALL ERRORS("FIELD_IO_LABEL_DERIVATIVE_INFO_GET",ERR,ERROR)
    CALL EXITS("FIELD_IO_LABEL_DERIVATIVE_INFO_GET")
  END FUNCTION FIELD_IO_LABEL_DERIVATIVE_INFO_GET    

  !
  !================================================================================================================================
  !  

  !>Get the field information 
  FUNCTION FIELD_IO_LABEL_FIELD_INFO_GET(COMPONENT, LABEL_TYPE, ERR, ERROR)
    !Argument variables   
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: COMPONENT
    INTEGER(INTG), INTENT(IN) :: LABEL_TYPE !<identitor for information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    TYPE(VARYING_STRING) :: FIELD_IO_LABEL_FIELD_INFO_GET
   
    CALL ENTERS("FIELD_IO_LABEL_FIELD_INFO_GET",ERR,ERROR,*999)    
   
    IF(.NOT.ASSOCIATED(COMPONENT)) THEN
       CALL FLAG_ERROR("component pointer in the input data is invalid",ERR,ERROR,*999)   
       GOTO 999        
    ENDIF
   
    IF(LABEL_TYPE/=FIELD_IO_FIELD_LABEL.AND.LABEL_TYPE/=FIELD_IO_VARIABLE_LABEL.AND.LABEL_TYPE/=FIELD_IO_COMPONENT_LABEL) THEN
       CALL FLAG_ERROR("indicator of type in input data is invalid",ERR,ERROR,*999)  
       GOTO 999         
    ENDIF
   
    FIELD=>COMPONENT%FIELD
    VARIABLE=>COMPONENT%FIELD_VARIABLE
    
    SELECT CASE(FIELD%TYPE)
      CASE(FIELD_GEOMETRIC_TYPE) !FIELD_GEOMETRIC_TYPE
        IF(LABEL_TYPE==FIELD_IO_FIELD_LABEL) THEN
           FIELD_IO_LABEL_FIELD_INFO_GET="field geometric type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
              !coordinate system              
              SELECT CASE (COMPONENT%REGION%COORDINATE_SYSTEM%TYPE)
                  CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
                    IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN                
                       FIELD_IO_LABEL_FIELD_INFO_GET="coordinates,  coordinate, rectangular cartesian"
                    ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                       IF(COMPONENT%COMPONENT_NUMBER==1) THEN
                          FIELD_IO_LABEL_FIELD_INFO_GET="x"
                       ELSE IF(COMPONENT%COMPONENT_NUMBER==2) THEN
                          FIELD_IO_LABEL_FIELD_INFO_GET="y"
                       ELSE IF(COMPONENT%COMPONENT_NUMBER==3) THEN
                          FIELD_IO_LABEL_FIELD_INFO_GET="z"
                       ENDIF                             
                    ENDIF                        
                  !CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
                  !CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
                  !CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
                  !CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
                  CASE DEFAULT
                   IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN                
                      FIELD_IO_LABEL_FIELD_INFO_GET="unknown" !coordinates, coordinate, rectangular cartesian,
                   ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                      FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
                   ENDIF                        
              END SELECT
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="Normal_derivative,  field,  normal derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="first_time_derivative,  field,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="second_time_derivative,  field,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown_geometry,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_IO_FIELD_LABEL   
      CASE(FIELD_FIBRE_TYPE)          
        IF(LABEL_TYPE==FIELD_IO_FIELD_LABEL) THEN
           FIELD_IO_LABEL_FIELD_INFO_GET="field fibres type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="fiber,  standand variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="norm_der_fiber,  normal derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="first_time_fiber,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="second_time_fiber,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown_fiber,  unknown field variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_IO_FIELD_LABEL           
      CASE(FIELD_GENERAL_TYPE)   
        IF(LABEL_TYPE==FIELD_IO_FIELD_LABEL) THEN
           FIELD_IO_LABEL_FIELD_INFO_GET="field general type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="general_variabe,  field,  string"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="norm_dev_variable,  field,  string"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="first_time_variable,  field,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="second_time_variable,  field,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown_general,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_IO_FIELD_LABEL           
      CASE(FIELD_MATERIAL_TYPE)      
        IF(LABEL_TYPE==FIELD_IO_FIELD_LABEL) THEN
           FIELD_IO_LABEL_FIELD_INFO_GET="field material type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="material,  field,  standand variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="normal_material,  field,  normal derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="fist_time_material,  field,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="second_time_material,  field,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown material,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_IO_FIELD_LABEL           
      CASE DEFAULT
        IF(LABEL_TYPE==FIELD_IO_FIELD_LABEL) THEN
           FIELD_IO_LABEL_FIELD_INFO_GET="unknown field type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown,  field,  unknown standand variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown,  field,  unknown normal derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown,  field,  unknown firt time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown, field,  unknown second time derivative of variable"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"Node:*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==FIELD_IO_VARIABLE_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET="unknown,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==FIELD_IO_COMPONENT_LABEL) THEN
                FIELD_IO_LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_IO_FIELD_LABEL           
    END SELECT
             
    CALL EXITS("FIELD_IO_LABEL_FIELD_INFO_GET")
    RETURN
999 CALL ERRORS("FIELD_IO_LABEL_FIELD_INFO_GET",ERR,ERROR)
    CALL EXITS("FIELD_IO_LABEL_FIELD_INFO_GET")
  END FUNCTION FIELD_IO_LABEL_FIELD_INFO_GET  

  !!
  !!================================================================================================================================
  !!  

  !!>Write the header of a group nodes using FORTRAIN  
  !SUBROUTINE FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN(LOCAL_PROCESS_NODAL_INFO_SET, LOCAL_NODAL_NUMBER, MAX_NUM_OF_NODAL_DERIVATIVES, &
  !&my_computational_node_number, FILE_ID, ERR,ERROR, *)
  !  !Argument variables  
  !  TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_NODAL_INFO_SET  !<LOCAL_PROCESS_NODAL_INFO_SET           
  !  INTEGER(INTG), INTENT(IN) :: LOCAL_NODAL_NUMBER !<LOCAL_NUMBER IN THE NODAL IO LIST
  !  INTEGER(INTG), INTENT(INOUT) :: MAX_NUM_OF_NODAL_DERIVATIVES !<MAX_NUM_OF_NODAL_DERIVATIVES
  !  INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<local process number    
  !  INTEGER(INTG), INTENT(IN) :: FILE_ID !< FILE ID    
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(FIELD_TYPE), POINTER :: field_ptr
  !  TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable_ptr   
  !  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
  !  TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
  !  TYPE(VARYING_STRING) :: LINE, LABEL
  !  INTEGER(INTG) :: NUM_OF_FIELDS, NUM_OF_VARIABLES, NUM_OF_NODAL_DEV
  !  INTEGER(INTG) :: domain_idx, domain_no, MY_DOMAIN_INDEX, local_number, global_number
  !  INTEGER(INTG), POINTER :: GROUP_FIELDS(:), GROUP_VARIABLES(:), GROUP_DERIVATIVES(:)
  !  INTEGER(INTG) :: field_idx, comp_idx, comp_idx1, value_idx, var_idx, global_var_idx !dev_idx,   
  ! 
  !  CALL ENTERS("FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN",ERR,ERROR,*999)    
  !  
  !  !colllect nodal header information for IO first
  !  
  !  !!get the number of this computational node from mpi pool
  !  !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
  !  !IF(ERR/=0) GOTO 999        
  !  
  !  !attach the temporary pointer
  !  !tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS    
  !  
  !  !collect maximum number of nodal derivatives, number of fields and variables 
  !  NUM_OF_FIELDS=0
  !  NUM_OF_VARIABLES=0
  !  MAX_NUM_OF_NODAL_DERIVATIVES=0
  !  global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_NODAL_NUMBER)
  !  NULLIFY(field_ptr)     
  !  NULLIFY(variable_ptr)
  !  DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
  !     !calculate the number of fields
  !     IF (.NOT.ASSOCIATED(field_ptr, target=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !          &COMPONENTS(comp_idx)%PTR%FIELD)) THEN
  !        NUM_OF_FIELDS=NUM_OF_FIELDS+1
  !        field_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS (comp_idx)%PTR%FIELD
  !     ENDIF         
  !     
  !     !calculate the number of variables
  !     IF (.NOT.ASSOCIATED(variable_ptr, target=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !          &COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE)) THEN
  !        NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
  !        variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE
  !     ENDIF      
  ! 
  !    !finding the local numbering through the global to local mapping   
  !     DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
  !         &COMPONENTS(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
  !     !get the domain index for this variable component according to my own computional node number
  !     DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
  !        domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
  !        IF(domain_no==my_computational_node_number) THEN
  !           MY_DOMAIN_INDEX=domain_idx
  !           EXIT !out of loop--domain_idx
  !        ENDIF
  !     ENDDO !domain_idX
  !     local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
  !     !use local domain information find the out the maximum number of derivatives
  !     DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES          
  !     MAX_NUM_OF_NODAL_DERIVATIVES=MAX(DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES,MAX_NUM_OF_NODAL_DERIVATIVES)                
  !  ENDDO !comp_idx         
  !  !Allocate the momery for group of field variables
  !  ALLOCATE(GROUP_FIELDS(NUM_OF_FIELDS),STAT=ERR)
  !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty field buffer in IO",ERR,ERROR,*999)
  !  !Allocate the momery for group of field components   
  !  ALLOCATE(GROUP_VARIABLES(NUM_OF_VARIABLES),STAT=ERR)
  !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
  !  !Allocate the momery for group of maximum number of derivatives     
  !  ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivatives buffer in IO",ERR,ERROR,*999)    
  !  
  !  !fill information into the group of fields and variables
  !  NUM_OF_FIELDS=0
  !  NUM_OF_VARIABLES=0   
  !  NULLIFY(field_ptr)     
  !  NULLIFY(variable_ptr)      
  !  GROUP_FIELDS(:)=0 !the item in this arrary is the number of variables in the same field
  !  GROUP_VARIABLES(:)=0 !the item in this arrary is the number of components in the same variable     
  !  DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
  !     !grouping field variables and components together
  !     IF((.NOT.ASSOCIATED(field_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !          &COMPONENTS(comp_idx)%PTR%FIELD)).AND.(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET% &
  !          &NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE))) THEN !different field and variables            
  !        !add one new variable
  !        NUM_OF_FIELDS=NUM_OF_FIELDS+1
  !        GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1 
  !        !add one new component
  !        NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
  !        GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
  !        field_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD
  !        variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE   
  !     ELSE IF (ASSOCIATED(field_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !        &COMPONENTS(comp_idx)%PTR%FIELD).AND.(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%&
  !        &NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE))) THEN !the same field and  different variables
  !        !add one new variable
  !        GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1
  !        !add one new component
  !        NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
  !        GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
  !        variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE   
  !     ELSE  !different components of the same variable
  !        !add one new component
  !        GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1                               
  !     ENDIF !field_ptr/=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS%COMPONENTS(comp_idx)%PTR%FIELD     
  !  ENDDO  !comp_idx        
  !  
  !  !write out the nodal header
  !  var_idx=1
  !  comp_idx=1
  !  field_idx=1
  !  value_idx=1
  !  comp_idx1=1
  !  global_var_idx=0
  !  
  !  CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
  !  IF(LINE/=" "//"#Fields="//TRIM(NUMBER_TO_VSTRING(SUM(GROUP_FIELDS(1:NUM_OF_FIELDS)),"*",ERR,ERROR))) &
  !     & CALL FLAG_ERROR("Fields number in the Header part do not match",ERR,ERROR,*999)
  !       
  !  DO field_idx=1, NUM_OF_FIELDS              
  !     !write out the field information
  !     !LABEL=FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx1)%PTR, FIELD_IO_FIELD_LABEL,ERR,ERROR)
  !     !IF(ERR/=0) THEN
  !     !   CALL FLAG_ERROR("can not get field label",ERR,ERROR,*999)     
  !     !ENDIF
  !     !CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)    
  !     !IF(LINE/=TRIM(NUMBER_TO_VSTRING(field_idx,"*",ERR,ERROR))//") "//TRIM(LABEL)&
  !     !&//" , #variables="//TRIM(NUMBER_TO_VSTRING(GROUP_FIELDS(field_idx),"*",ERR,ERROR)) &       
  !     !CALL FLAG_ERROR("Variable number in the Header part do not match",ERR,ERROR,*999)              
  !     
  !     DO var_idx=1, GROUP_FIELDS(field_idx)
  !        global_var_idx=global_var_idx+1
  !        !write out the field information
  !        LABEL="  "//TRIM(NUMBER_TO_VSTRING(global_var_idx,"*",ERR,ERROR))//") "&
  !        &//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
  !        &COMPONENTS(comp_idx1)%PTR, FIELD_IO_VARIABLE_LABEL,ERR,ERROR)
  !        IF(ERR/=0) THEN
  !           CALL FLAG_ERROR("can not get variable label",ERR,ERROR,*999)     
  !           GOTO 999               
  !        ENDIF        
  !        CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
  !        IF(LINE/=TRIM(LABEL)//", #Components="//TRIM(NUMBER_TO_VSTRING(GROUP_VARIABLES(global_var_idx),"*",ERR,ERROR))) &
  !           & CALL FLAG_ERROR("Components number in the Header part do not match",ERR,ERROR,*999)       
  !                      
  !        
  !        DO comp_idx=1, GROUP_VARIABLES(global_var_idx)
  !           !write out the component information
  !           LABEL="   "//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
  !           &COMPONENTS(comp_idx1)%PTR, FIELD_IO_COMPONENT_LABEL,ERR,ERROR)
  !           IF(ERR/=0) THEN
  !              CALL FLAG_ERROR("can not get component label",ERR,ERROR,*999)     
  !              GOTO 999               
  !           ENDIF        
  !           LINE=TRIM(LABEL)//"."                          
  !           
  !           !finding the local numbering through the global to local mapping   
  !           DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%&
  !              &DOMAIN%MAPPINGS%NODES       
  !           !get the domain index for this variable component according to my own computional node number
  !           DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
  !              domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
  !              IF(domain_no==my_computational_node_number) THEN
  !                 MY_DOMAIN_INDEX=domain_idx
  !                 EXIT !out of loop--domain_idx
  !              ENDIF
  !           ENDDO !domain_idX
  !           local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
  !           !use local domain information find the out the maximum number of derivatives
  !           DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
  !           !get the nodal partial derivatives
  !           NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
  !           GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
  !           !sort  the partial derivatives
  !           CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
  !           !get the derivative name             
  !           LABEL=FIELD_IO_LABEL_DERIVATIVE_INFO_GET(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV), NUM_OF_NODAL_DEV, &
  !              &FIELD_IO_DERIVATIVE_LABEL,ERR,ERROR)
  !           IF(ERR/=0) THEN
  !              CALL FLAG_ERROR("can not get derivative label",ERR,ERROR,*999)     
  !              GOTO 999               
  !           ENDIF
  !           !write out the header             
  !           CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
  !           !assemble the header        
  !           IF(LINE/=LINE//"  Value index= "//TRIM(NUMBER_TO_VSTRING(value_idx,"*",ERR,ERROR))&
  !              &//", #Derivatives= "//TRIM(NUMBER_TO_VSTRING(NUM_OF_NODAL_DEV-1,"*",ERR,ERROR))//TRIM(LABEL)) &
  !              & CALL FLAG_ERROR("Value index in the Header part do not match",ERR,ERROR,*999)
  !           !increase the component index              
  !           comp_idx1=comp_idx1+1
  !           !increase the value index
  !           value_idx=value_idx+NUM_OF_NODAL_DEV
  !        ENDDO !comp_idx
  !     ENDDO !var_idx
  !  ENDDO !field_idx
  !
  !  !release temporary memory
  !  IF(ASSOCIATED(GROUP_FIELDS)) DEALLOCATE(GROUP_FIELDS)
  !  IF(ASSOCIATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)
  !  IF(ASSOCIATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)    
  !      
  !  CALL EXITS("FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN")
  !  RETURN 1       
  !END SUBROUTINE FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN  

  !
  !================================================================================================================================
  !  

  !>Write the header of a group nodes using FORTRAIN  
  SUBROUTINE FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN(LOCAL_PROCESS_NODAL_INFO_SET, LOCAL_NODAL_NUMBER, MAX_NUM_OF_NODAL_DERIVATIVES, &
  &my_computational_node_number, FILE_ID, ERR,ERROR, *)
    !Argument variables  
    TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_NODAL_INFO_SET  !<LOCAL_PROCESS_NODAL_INFO_SET           
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODAL_NUMBER !<LOCAL_NUMBER IN THE NODAL IO LIST
    INTEGER(INTG), INTENT(INOUT) :: MAX_NUM_OF_NODAL_DERIVATIVES !<MAX_NUM_OF_NODAL_DERIVATIVES
    INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN) :: FILE_ID !< FILE ID    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: field_ptr
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable_ptr   
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
    TYPE(VARYING_STRING) :: LINE, LABEL
    INTEGER(INTG), ALLOCATABLE :: GROUP_FIELDS(:), GROUP_VARIABLES(:), GROUP_DERIVATIVES(:)    
    INTEGER(INTG) :: NUM_OF_FIELDS, NUM_OF_VARIABLES, NUM_OF_NODAL_DEV
    INTEGER(INTG) :: domain_idx, domain_no, MY_DOMAIN_INDEX, local_number, global_number    
    INTEGER(INTG) :: field_idx, comp_idx, comp_idx1, value_idx, var_idx, global_var_idx !dev_idx,   
   
    CALL ENTERS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN",ERR,ERROR,*999)    
    
    !colllect nodal header information for IO first
    
    !!get the number of this computational node from mpi pool
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
    !IF(ERR/=0) GOTO 999        
    
    !attach the temporary pointer
    !tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS    
    
    !collect maximum number of nodal derivatives, number of fields and variables 
    NUM_OF_FIELDS=0
    NUM_OF_VARIABLES=0
    MAX_NUM_OF_NODAL_DERIVATIVES=0
    global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_NODAL_NUMBER)
    NULLIFY(field_ptr)     
    NULLIFY(variable_ptr)
    DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of fields
       IF (.NOT.ASSOCIATED(field_ptr, target=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
            &COMPONENTS(comp_idx)%PTR%FIELD)) THEN
          NUM_OF_FIELDS=NUM_OF_FIELDS+1
          field_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS (comp_idx)%PTR%FIELD
       ENDIF         
       
       !calculate the number of variables
       IF (.NOT.ASSOCIATED(variable_ptr, target=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
            &COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE)) THEN
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE
       ENDIF      
   
      !finding the local numbering through the global to local mapping   
       DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
           &COMPONENTS(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
       !get the domain index for this variable component according to my own computional node number
       DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          IF(domain_no==my_computational_node_number) THEN
             MY_DOMAIN_INDEX=domain_idx
             EXIT !out of loop--domain_idx
          ENDIF
       ENDDO !domain_idX
       local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
       !use local domain information find the out the maximum number of derivatives
       DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES          
       MAX_NUM_OF_NODAL_DERIVATIVES=MAX(DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES,MAX_NUM_OF_NODAL_DERIVATIVES)                
    ENDDO !comp_idx         
    !Allocate the momery for group of field variables
    ALLOCATE(GROUP_FIELDS(NUM_OF_FIELDS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty field buffer in IO",ERR,ERROR,*999)
    !Allocate the momery for group of field components   
    ALLOCATE(GROUP_VARIABLES(NUM_OF_VARIABLES),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
    !Allocate the momery for group of maximum number of derivatives     
    ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivatives buffer in IO",ERR,ERROR,*999)    
    
    !fill information into the group of fields and variables
    NUM_OF_FIELDS=0
    NUM_OF_VARIABLES=0   
    NULLIFY(field_ptr)     
    NULLIFY(variable_ptr)      
    GROUP_FIELDS(:)=0 !the item in this arrary is the number of variables in the same field
    GROUP_VARIABLES(:)=0 !the item in this arrary is the number of components in the same variable     
    DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !grouping field variables and components together
       IF((.NOT.ASSOCIATED(field_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
            &COMPONENTS(comp_idx)%PTR%FIELD)).AND.(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET% &
            &NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE))) THEN !different field and variables            
          !add one new variable
          NUM_OF_FIELDS=NUM_OF_FIELDS+1
          GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1 
          !add one new component
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
          field_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD
          variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE   
       ELSE IF (ASSOCIATED(field_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
          &COMPONENTS(comp_idx)%PTR%FIELD).AND.(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%&
          &NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE))) THEN !the same field and  different variables
          !add one new variable
          GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1
          !add one new component
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
          variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE   
       ELSE  !different components of the same variable
          !add one new component
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1                               
       ENDIF !field_ptr/=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS%COMPONENTS(comp_idx)%PTR%FIELD     
    ENDDO  !comp_idx        
    
    !write out the nodal header
    var_idx=1
    comp_idx=1
    field_idx=1
    value_idx=1
    comp_idx1=1
    global_var_idx=0
    LINE=" #Fields="//TRIM(NUMBER_TO_VSTRING(SUM(GROUP_FIELDS(1:NUM_OF_FIELDS)),"*",ERR,ERROR))
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    DO field_idx=1, NUM_OF_FIELDS              
       !write out the field information
       !LABEL=FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx1)%PTR, FIELD_IO_FIELD_LABEL,ERR,ERROR)
       !IF(ERR/=0) THEN
       !   CALL FLAG_ERROR("can not get field label",ERR,ERROR,*999)     
       !   GOTO 999               
       !ENDIF        
       !LINE=TRIM(NUMBER_TO_VSTRING(field_idx,"*",ERR,ERROR))//") "//TRIM(LABEL)&
       !&//" , #variables="//TRIM(NUMBER_TO_VSTRING(GROUP_FIELDS(field_idx),"*",ERR,ERROR))       
       !CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
       
       DO var_idx=1, GROUP_FIELDS(field_idx)
          global_var_idx=global_var_idx+1
          !write out the field information
          LABEL=" "//TRIM(NUMBER_TO_VSTRING(global_var_idx,"*",ERR,ERROR))//") "&
          &//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
          &COMPONENTS(comp_idx1)%PTR, FIELD_IO_VARIABLE_LABEL,ERR,ERROR)
          IF(ERR/=0) THEN
             CALL FLAG_ERROR("can not get variable label",ERR,ERROR,*999)     
             GOTO 999               
          ENDIF        
          LINE=TRIM(LABEL)//", #Components="//TRIM(NUMBER_TO_VSTRING(GROUP_VARIABLES(global_var_idx),"*",ERR,ERROR))        
          CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
          
          DO comp_idx=1, GROUP_VARIABLES(global_var_idx)
             !write out the component information
             LABEL="   "//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
             &COMPONENTS(comp_idx1)%PTR, FIELD_IO_COMPONENT_LABEL,ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get component label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF        
             LINE=TRIM(LABEL)//"."                          
             
             !finding the local numbering through the global to local mapping   
             DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%&
                &DOMAIN%MAPPINGS%NODES       
             !get the domain index for this variable component according to my own computional node number
             DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
                domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
                IF(domain_no==my_computational_node_number) THEN
                   MY_DOMAIN_INDEX=domain_idx
                   EXIT !out of loop--domain_idx
                ENDIF
             ENDDO !domain_idX
             local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
             !use local domain information find the out the maximum number of derivatives
             DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
             !get the nodal partial derivatives
             NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
             GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
             !sort  the partial derivatives
             CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
             !get the derivative name             
             LABEL=FIELD_IO_LABEL_DERIVATIVE_INFO_GET(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV), NUM_OF_NODAL_DEV, &
                &FIELD_IO_DERIVATIVE_LABEL,ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get derivative label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF
             !assemble the header        
             LINE=LINE//"  Value index= "//TRIM(NUMBER_TO_VSTRING(value_idx,"*",ERR,ERROR))&
             &//", #Derivatives= "//TRIM(NUMBER_TO_VSTRING(NUM_OF_NODAL_DEV-1,"*",ERR,ERROR))//TRIM(LABEL)
             !write out the header             
             CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)
             !increase the component index              
             comp_idx1=comp_idx1+1
             !increase the value index
             value_idx=value_idx+NUM_OF_NODAL_DEV
          ENDDO !comp_idx
       ENDDO !var_idx
    ENDDO !field_idx

    !release temporary memory
    IF(ALLOCATED(GROUP_FIELDS)) DEALLOCATE(GROUP_FIELDS)
    IF(ALLOCATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)
    IF(ALLOCATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)    
        
    CALL EXITS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN")
    RETURN
999 CALL ERRORS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN",ERR,ERROR)
    CALL EXITS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN")
    RETURN 1       
  END SUBROUTINE FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN  

  !!
  !!================================================================================================================================
  !!  
  !
  !!>Write the header of a group nodes using MPIIO  
  !SUBROUTINE FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO(LOCAL_PROCESS_NODAL_INFO_SET, LOCAL_NODAL_NUMBER, MAX_NUM_OF_NODAL_DERIVATIVES, &
  !&my_computational_node_number, FILE_ID, ERR,ERROR, *)
  !  !Argument variables  
  !  TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT) :: LOCAL_PROCESS_NODAL_INFO_SET  !<LOCAL_PROCESS_NODAL_INFO_SET           
  !  INTEGER(INTG), INTENT(IN) :: LOCAL_NODAL_NUMBER !<LOCAL_NUMBER IN THE NODAL IO LIST
  !  INTEGER(INTG), INTENT(INOUT) :: MAX_NUM_OF_NODAL_DERIVATIVES !<MAX_NUM_OF_NODAL_DERIVATIVES
  !  INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<local process number    
  !  INTEGER(INTG), INTENT(IN) :: FILE_ID !< FILE ID    
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(FIELD_TYPE), POINTER :: field_ptr
  !  TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable_ptr   
  !  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
  !  TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
  !  TYPE(VARYING_STRING) :: LINE, LABEL
  !  INTEGER(INTG) :: NUM_OF_FIELDS, NUM_OF_VARIABLES, NUM_OF_NODAL_DEV
  !  INTEGER(INTG) :: domain_idx, domain_no, MY_DOMAIN_INDEX, local_number, global_number
  !  INTEGER(INTG), POINTER :: GROUP_FIELDS(:), GROUP_VARIABLES(:), GROUP_DERIVATIVES(:)
  !  INTEGER(INTG) :: field_idx, comp_idx, comp_idx1, value_idx, var_idx, global_var_idx !dev_idx,   
  ! 
  !  CALL ENTERS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO",ERR,ERROR,*999)    
  !  
  !  !colllect nodal header information for IO first
  !  
  !  !!get the number of this computational node from mpi pool
  !  !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
  !  !IF(ERR/=0) GOTO 999        
  !  
  !  !attach the temporary pointer
  !  !tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS    
  !  
  !  !collect maximum number of nodal derivatives, number of fields and variables 
  !  NUM_OF_FIELDS=0
  !  NUM_OF_VARIABLES=0
  !  MAX_NUM_OF_NODAL_DERIVATIVES=0
  !  global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_NODAL_NUMBER)
  !  NULLIFY(field_ptr)     
  !  NULLIFY(variable_ptr)
  !  DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
  !     !calculate the number of fields
  !     IF (.NOT.ASSOCIATED(field_ptr, target=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !          &COMPONENTS(comp_idx)%PTR%FIELD)) THEN
  !        NUM_OF_FIELDS=NUM_OF_FIELDS+1
  !        field_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS (comp_idx)%PTR%FIELD
  !     ENDIF         
  !     
  !     !calculate the number of variables
  !     IF (.NOT.ASSOCIATED(variable_ptr, target=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !          &COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE)) THEN
  !        NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
  !        variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE
  !     ENDIF      
  ! 
  !    !finding the local numbering through the global to local mapping   
  !     DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
  !         &COMPONENTS(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
  !     !get the domain index for this variable component according to my own computional node number
  !     DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
  !        domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
  !        IF(domain_no==my_computational_node_number) THEN
  !           MY_DOMAIN_INDEX=domain_idx
  !           EXIT !out of loop--domain_idx
  !        ENDIF
  !     ENDDO !domain_idX
  !     local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
  !     !use local domain information find the out the maximum number of derivatives
  !     DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES          
  !     MAX_NUM_OF_NODAL_DERIVATIVES=MAX(DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES,MAX_NUM_OF_NODAL_DERIVATIVES)                
  !  ENDDO !comp_idx         
  !  !Allocate the momery for group of field variables
  !  ALLOCATE(GROUP_FIELDS(NUM_OF_FIELDS),STAT=ERR)
  !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty field buffer in IO",ERR,ERROR,*999)
  !  !Allocate the momery for group of field components   
  !  ALLOCATE(GROUP_VARIABLES(NUM_OF_VARIABLES),STAT=ERR)
  !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
  !  !Allocate the momery for group of maximum number of derivatives     
  !  ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivatives buffer in IO",ERR,ERROR,*999)    
  !  
  !  !fill information into the group of fields and variables
  !  NUM_OF_FIELDS=0
  !  NUM_OF_VARIABLES=0   
  !  NULLIFY(field_ptr)     
  !  NULLIFY(variable_ptr)      
  !  GROUP_FIELDS(:)=0 !the item in this arrary is the number of variables in the same field
  !  GROUP_VARIABLES(:)=0 !the item in this arrary is the number of components in the same variable     
  !  DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
  !     !grouping field variables and components together
  !     IF((.NOT.ASSOCIATED(field_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !          &COMPONENTS(comp_idx)%PTR%FIELD)).AND.(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET% &
  !          &NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE))) THEN !different field and variables            
  !        !add one new variable
  !        NUM_OF_FIELDS=NUM_OF_FIELDS+1
  !        GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1 
  !        !add one new component
  !        NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
  !        GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
  !        field_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD
  !        variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE   
  !     ELSE IF (ASSOCIATED(field_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)% &
  !        &COMPONENTS(comp_idx)%PTR%FIELD).AND.(.NOT.ASSOCIATED(variable_ptr,TARGET=LOCAL_PROCESS_NODAL_INFO_SET%&
  !        &NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE))) THEN !the same field and  different variables
  !        !add one new variable
  !        GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1
  !        !add one new component
  !        NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
  !        GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
  !        variable_ptr=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%FIELD_VARIABLE   
  !     ELSE  !different components of the same variable
  !        !add one new component
  !        GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1                               
  !     ENDIF !field_ptr/=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS%COMPONENTS(comp_idx)%PTR%FIELD     
  !  ENDDO  !comp_idx        
  !  
  !  !write out the nodal header
  !  var_idx=1
  !  comp_idx=1
  !  field_idx=1
  !  value_idx=1
  !  comp_idx1=1
  !  global_var_idx=0
  !  LINE=" "//"#Fields="//TRIM(NUMBER_TO_VSTRING(SUM(GROUP_FIELDS(1:NUM_OF_FIELDS)),"*",ERR,ERROR))
  !  CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
  !  DO field_idx=1, NUM_OF_FIELDS              
  !     !write out the field information
  !     !LABEL=FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx1)%PTR, FIELD_IO_FIELD_LABEL,ERR,ERROR)
  !     !IF(ERR/=0) THEN
  !     !   CALL FLAG_ERROR("can not get field label",ERR,ERROR,*999)     
  !     !   GOTO 999               
  !     !ENDIF        
  !     !LINE=TRIM(NUMBER_TO_VSTRING(field_idx,"*",ERR,ERROR))//") "//TRIM(LABEL)&
  !     !&//" , #variables="//TRIM(NUMBER_TO_VSTRING(GROUP_FIELDS(field_idx),"*",ERR,ERROR))       
  !     !CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
  !     
  !     DO var_idx=1, GROUP_FIELDS(field_idx)
  !        global_var_idx=global_var_idx+1
  !        !write out the field information
  !        LABEL="  "//TRIM(NUMBER_TO_VSTRING(global_var_idx,"*",ERR,ERROR))//") "&
  !        &//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
  !        &COMPONENTS(comp_idx1)%PTR, FIELD_IO_VARIABLE_LABEL,ERR,ERROR)
  !        IF(ERR/=0) THEN
  !           CALL FLAG_ERROR("can not get variable label",ERR,ERROR,*999)     
  !           GOTO 999               
  !        ENDIF        
  !        LINE=TRIM(LABEL)//", #Components="//TRIM(NUMBER_TO_VSTRING(GROUP_VARIABLES(global_var_idx),"*",ERR,ERROR))        
  !        CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
  !        
  !        DO comp_idx=1, GROUP_VARIABLES(global_var_idx)
  !           !write out the component information
  !           LABEL="   "//FIELD_IO_LABEL_FIELD_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%&
  !           &COMPONENTS(comp_idx1)%PTR, FIELD_IO_COMPONENT_LABEL,ERR,ERROR)
  !           IF(ERR/=0) THEN
  !              CALL FLAG_ERROR("can not get component label",ERR,ERROR,*999)     
  !              GOTO 999               
  !           ENDIF        
  !           LINE=TRIM(LABEL)//"."                          
  !           
  !           !finding the local numbering through the global to local mapping   
  !           DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%&
  !              &DOMAIN%MAPPINGS%NODES       
  !           !get the domain index for this variable component according to my own computional node number
  !           DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
  !              domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
  !              IF(domain_no==my_computational_node_number) THEN
  !                 MY_DOMAIN_INDEX=domain_idx
  !                 EXIT !out of loop--domain_idx
  !              ENDIF
  !           ENDDO !domain_idX
  !           local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
  !           !use local domain information find the out the maximum number of derivatives
  !           DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
  !           !get the nodal partial derivatives
  !           NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
  !           GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
  !           !sort  the partial derivatives
  !           CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
  !           !get the derivative name             
  !           LABEL=FIELD_IO_LABEL_DERIVATIVE_INFO_GET(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV), NUM_OF_NODAL_DEV, &
  !              &FIELD_IO_DERIVATIVE_LABEL,ERR,ERROR)
  !           IF(ERR/=0) THEN
  !              CALL FLAG_ERROR("can not get derivative label",ERR,ERROR,*999)     
  !              GOTO 999               
  !           ENDIF
  !           !assemble the header        
  !           LINE=LINE//"  Value index= "//TRIM(NUMBER_TO_VSTRING(value_idx,"*",ERR,ERROR))&
  !           &//", #Derivatives= "//TRIM(NUMBER_TO_VSTRING(NUM_OF_NODAL_DEV-1,"*",ERR,ERROR))//TRIM(LABEL)
  !           !write out the header             
  !           CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)
  !           !increase the component index              
  !           comp_idx1=comp_idx1+1
  !           !increase the value index
  !           value_idx=value_idx+NUM_OF_NODAL_DEV
  !        ENDDO !comp_idx
  !     ENDDO !var_idx
  !  ENDDO !field_idx!
  !
  !  !release temporary memory
  !  IF(ASSOCIATED(GROUP_FIELDS)) DEALLOCATE(GROUP_FIELDS)
  !  IF(ASSOCIATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)
  !  IF(ASSOCIATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)    
  !      
  !  CALL EXITS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO")
  !  RETURN 1       
  !END SUBROUTINE FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO  

  !!
  !!================================================================================================================================
  !!  
  !
  !!>Write all the nodal information from LOCAL_PROCESS_NODAL_INFO_SET to single exnode files       
  !SUBROUTINE FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE(LOCAL_PROCESS_NODAL_INFO_SET, NAME, my_computational_node_number, &
  !  &computational_node_numbers,ERR, ERROR, *)
  !  !the reason that my_computational_node_number is used in the argument is for future extension
  !  !Argument variables   
  !  TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT):: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
  !  TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
  !  INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
  !  INTEGER(INTG), INTENT(IN):: computational_node_numbers !<total process number       
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(VARYING_STRING) :: LINE 
  !  TYPE(VARYING_STRING) :: FILE_NAME !the prefix name of file.    
  !  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
  !  TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
  !  INTEGER(INTG) :: bufsize
  !  INTEGER(INTG) :: FILE_ID, domain_idx, domain_no, local_number, global_number, MY_DOMAIN_INDEX
  !  INTEGER(INTG), ALLOCATABLE :: GROUP_DERIVATIVES(:)    
  !  INTEGER(INTG) :: nn, comp_idx,  dev_idx, NUM_OF_NODAL_DEV, MAX_NUM_OF_NODAL_DERIVATIVES
  !  REAL(DP), ALLOCATABLE :: NODAL_BUFFER(:)
  !  REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
  !
  !  CALL ENTERS("FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE",ERR,ERROR,*999)    
  !  
  !  !get my own computianal node number--be careful the rank of process in the MPI pool 
  !  !is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
  !  !will be a secured way to get the number         
  !  FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(my_computational_node_number,"*",ERR,ERROR))//".exnode"        
  !  FILE_ID=1029
  !  MAX_NUM_OF_NODAL_DERIVATIVES=0
  !  bufsize=999999
  !  
  !  IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
  !     CALL FLAG_ERROR("the nodal information set in input is invalid",ERR,ERROR,*999)            
  !  ENDIF
  !
  !  IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
  !     CALL FLAG_ERROR("the nodal information set is not associated with any numbering list",ERR,ERROR,*999)            
  !  ENDIF
  !
  !  IF(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0) THEN 
  !     CALL FLAG_ERROR("the nodal information set does not contain any nodes",ERR,ERROR,*999)            
  !  ENDIF
  !
  !  IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(1)%SAME_HEADER==.TRUE.) THEN
  !     CALL FLAG_ERROR("the first header flag of nodal information set should be false",ERR,ERROR,*999)            
  !  ENDIF
  !  
  !  !NULLIFY(NODAL_BUFFER)    
  !  !NULLIFY(GROUP_DERIVATIVES)
  !  
  !  !open a file 
  !  CALL FIELD_IO_MPIIO_FILE_OPEN(FILE_ID, FILE_NAME, my_computational_node_number, bufsize, ERR,ERROR,*999)
  !  
  !  !write out the group name    
  !  LINE=FIELD_IO_FILEDS_GROUP_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS, ERR,ERROR)    
  !  IF(ERR/=0) THEN
  !     CALL FLAG_ERROR("can not get group namein IO",ERR,ERROR,*999)     
  !     GOTO 999               
  !  ENDIF               
  !  CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
  !  !write out the number of files    
  !  LINE=FIELD_IO_MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)      
  !  IF(ERR/=0) THEN
  !     CALL FLAG_ERROR("can not get multiple file information in IO",ERR,ERROR,*999)     
  !     GOTO 999               
  !  ENDIF             
  !  !CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
  !
  !  DO nn=1, LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
  !    
  !     !tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS         
  !     global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)
  !     
  !     !check whether need to write out the nodal information header  
  !     IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.) THEN
  !      !write out the nodal header
  !         
  !        CALL FIELD_IO_EXPORT_NODAL_GROUP_HEADER_MPIIO(LOCAL_PROCESS_NODAL_INFO_SET, nn, MAX_NUM_OF_NODAL_DERIVATIVES, &
  !        &my_computational_node_number, FILE_ID, ERR,ERROR,*999) 
  !         !value_idx=value_idx-1 !the len of NODAL_BUFFER                      
  !        !checking: whether need to allocate temporary memory for Io writing
  !        IF(ALLOCATED(NODAL_BUFFER)) THEN
  !           IF(SIZE(NODAL_BUFFER)<MAX_NUM_OF_NODAL_DERIVATIVES) THEN 
  !              DEALLOCATE(NODAL_BUFFER)
  !              DEALLOCATE(GROUP_DERIVATIVES)
  !              ALLOCATE(NODAL_BUFFER(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty nodal buffer in IO writing",ERR,ERROR,*999)
  !              ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
  !           ENDIF
  !        ELSE
  !           ALLOCATE(NODAL_BUFFER(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty nodal buffer in IO writing",ERR,ERROR,*999)
  !           ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
  !        ENDIF
  !     ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.
  !                                 
  !     !write out the components' values of this node in this domain
  !     DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS   
  !
  !       !finding the local numbering through the global to local mapping   
  !        DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
  !        !get the domain index for this variable component according to my own computional node number
  !        DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
  !           domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
  !           IF(domain_no==my_computational_node_number) THEN
  !              MY_DOMAIN_INDEX=domain_idx
  !              EXIT !out of loop--domain_idx
  !           ENDIF
  !        ENDDO !domain_idX
  !        local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
  !        !use local domain information find the out the maximum number of derivatives
  !        DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
  !        
  !        !write out the user numbering of node if comp_idx ==1
  !        IF(comp_idx==1) THEN 
  !           LINE="Node:     "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NODES%NODES(local_number)%USER_NUMBER,"*",ERR,ERROR))
  !           CALL FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
  !        ENDIF
  !        
  !        !get the nodal partial derivatives
  !        NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
  !        GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
  !        !sort  the partial derivatives
  !        CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
  !        DO dev_idx=1, NUM_OF_NODAL_DEV
  !           NULLIFY(GEOMETRIC_PARAMETERS)
  !           IF(comp_idx==1) THEN
  !              CALL FIELD_PARAMETER_SET_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%FIELD,&
  !                                        &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
  !           ELSE IF (ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN, &
  !              &TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS (comp_idx-1)%PTR%DOMAIN)) THEN
  !              CALL FIELD_PARAMETER_SET_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%FIELD,&
  !                                        &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
  !           ENDIF                                       
  !           NODAL_BUFFER(dev_idx)=GEOMETRIC_PARAMETERS(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%&
  !             &PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(GROUP_DERIVATIVES(dev_idx),local_number,0))
  !           CALL FIELD_IO_MPIIO_FILE_WRITE_DP(FILE_ID, NODAL_BUFFER, NUM_OF_NODAL_DEV, ERR,ERROR,*999)                
  !        ENDDO
  !     ENDDO !comp_idx                   
  !     
  !  ENDDO!nn    
  !  
  !  !close a file    
  !  CALL FIELD_IO_MPIIO_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
  !
  !  !release the temporary memory
  !  IF(ALLOCATED(NODAL_BUFFER)) THEN
  !     DEALLOCATE(NODAL_BUFFER)
  !  ENDIF
  !  IF(ALLOCATED(GROUP_DERIVATIVES)) THEN
  !     DEALLOCATE(GROUP_DERIVATIVES)
  !  ENDIF
  !  
  !  CALL EXITS("FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE")
  !  RETURN 1  
  !END SUBROUTINE FIELD_IO_EXPORT_NODES_INTO_SINGLE_FILE
 
  !!
  !!================================================================================================================================
  !!  
  !
  !!>Read all the nodal information from LOCAL_PROCESS_NODAL_INFO_SET from local exnode files       
  !SUBROUTINE FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE(LOCAL_PROCESS_NODAL_INFO_SET, NAME, my_computational_node_number, &
  !  &computational_node_numbers,ERR, ERROR, *)
  !  !the reason that my_computational_node_number is used in the argument is for future extension
  !  !Argument variables   
  !  TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT):: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
  !  TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
  !  INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
  !  INTEGER(INTG), INTENT(IN):: computational_node_numbers !<total process number       
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(VARYING_STRING) :: LINE 
  !  TYPE(VARYING_STRING) :: FILE_NAME !the prefix name of file.    
  !  TYPE(VARYING_STRING) :: FILE_STATUS !the status for opening file.    
  !  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
  !  TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
  !  INTEGER(INTG) :: FILE_ID, domain_idx, domain_no, local_number, global_number, MY_DOMAIN_INDEX
  !  INTEGER(INTG), ALLOCATABLE :: GROUP_DERIVATIVES(:)    
  !  INTEGER(INTG) :: nn, comp_idx,  dev_idx, NUM_OF_NODAL_DEV, MAX_NUM_OF_NODAL_DERIVATIVES
  !  REAL(DP), ALLOCATABLE :: NODAL_BUFFER(:)
  !  REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
  !
  !  CALL ENTERS("FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE",ERR,ERROR,*999)    
  !  
  !  !get my own computianal node number--be careful the rank of process in the MPI pool 
  !  !is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
  !  !will be a secured way to get the number         
  !  FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(my_computational_node_number,"*",ERR,ERROR))//".exnode"        
  !  FILE_ID=1029
  !  MAX_NUM_OF_NODAL_DERIVATIVES=0
  !  
  !  IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
  !     CALL FLAG_ERROR("the nodal information set in input is invalid",ERR,ERROR,*999)            
  !  ENDIF
  !
  !  IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
  !     CALL FLAG_ERROR("the nodal information set is not associated with any numbering list",ERR,ERROR,*999)            
  !  ENDIF
  !
  !  IF(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0) THEN 
  !     CALL FLAG_ERROR("the nodal information set does not contain any nodes",ERR,ERROR,*999)            
  !  ENDIF
  !
  !  IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(1)%SAME_HEADER==.TRUE.) THEN
  !     CALL FLAG_ERROR("the first header flag of nodal information set should be false",ERR,ERROR,*999)            
  !  ENDIF
  !  
  !  !NULLIFY(NODAL_BUFFER)    
  !  !NULLIFY(GROUP_DERIVATIVES)
  !  
  !  !open a file for reading
  !  FILE_STATUS="OLD"
  !  CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
  !  IF(ERR/=0) THEN
  !     CALL FLAG_ERROR("can not find the file",ERR,ERROR,*999)     
  !     GOTO 999               
  !  ENDIF             
  !        
  !  !read out the group name    
  !  CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
  !  IF(LINE/=FIELD_IO_FILEDS_GROUP_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS, ERR,ERROR)) THEN   
  !     CALL FLAG_ERROR("group name does not match",ERR,ERROR,*999)     
  !     GOTO 999               
  !  ENDIF               
  !         
  !  !read out the number of files    
  !  LINE=FIELD_IO_MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)      
  !  IF(ERR/=0) THEN
  !     CALL FLAG_ERROR("can not get multiple file information in IO",ERR,ERROR,*999)     
  !     GOTO 999               
  !  ENDIF             
  !  !CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)       
  !
  !  DO nn=1, LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
  !    
  !     !tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS         
  !     global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)
  !     
  !     !check whether need to write out the nodal information header  
  !     IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.) THEN
  !      !write out the nodal header
  !         
  !        CALL FIELD_IO_IMPORT_NODAL_GROUP_HEADER_FORTRAN(LOCAL_PROCESS_NODAL_INFO_SET, nn, MAX_NUM_OF_NODAL_DERIVATIVES, &
  !        &my_computational_node_number, FILE_ID, ERR,ERROR,*999) 
  !         !value_idx=value_idx-1 !the len of NODAL_BUFFER                      
  !        !checking: whether need to allocate temporary memory for Io writing
  !        IF(ALLOCATED(NODAL_BUFFER)) THEN
  !           IF(SIZE(NODAL_BUFFER)<MAX_NUM_OF_NODAL_DERIVATIVES) THEN 
  !              DEALLOCATE(NODAL_BUFFER)
  !              DEALLOCATE(GROUP_DERIVATIVES)
  !              ALLOCATE(NODAL_BUFFER(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty nodal buffer in IO writing",ERR,ERROR,*999)
  !              ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
  !           ENDIF
  !        ELSE
  !           ALLOCATE(NODAL_BUFFER(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty nodal buffer in IO writing",ERR,ERROR,*999)
  !           ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
  !           IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
  !        ENDIF
  !     ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.
  !                                 
  !     !write out the components' values of this node in this domain
  !     DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS   
  !
  !       !finding the local numbering through the global to local mapping   
  !        DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
  !        !get the domain index for this variable component according to my own computional node number
  !        DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
  !           domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
  !           IF(domain_no==my_computational_node_number) THEN
  !              MY_DOMAIN_INDEX=domain_idx
  !              EXIT !out of loop--domain_idx
  !           ENDIF
  !        ENDDO !domain_idX
  !        local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
  !        !use local domain information find the out the maximum number of derivatives
  !        DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
  !        
  !        !write out the user numbering of node if comp_idx ==1
  !        IF(comp_idx==1) THEN 
  !           CALL FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, LINE, FILE_END, ERR,ERROR, *999)
  !           IF(LINE/="Node:     "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NODES%NODES(local_number)%USER_NUMBER,"*",ERR,ERROR))) THEN
  !              CALL FLAG_ERROR("node numbering does not match",ERR,ERROR,*999)
  !              GOTO 999      
  !           ENDIF                                      
  !        ENDIF
  !        
  !        !get the nodal partial derivatives
  !        NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
  !        GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
  !        !sort  the partial derivatives
  !        CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
  !        DO dev_idx=1, NUM_OF_NODAL_DEV
  !           CALL FIELD_IO_FORTRAN_FILE_READ_DP(FILE_ID, NODAL_BUFFER, NUM_OF_NODAL_DEV, ERR,ERROR,*999)              
  !           NULLIFY(GEOMETRIC_PARAMETERS)
  !           IF(comp_idx==1) THEN
  !              CALL FIELD_PARAMETER_SET_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%FIELD,&
  !                                        &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
  !           ELSE IF (ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN, &
  !              &TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS (comp_idx-1)%PTR%DOMAIN)) THEN
  !              CALL FIELD_PARAMETER_SET_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%FIELD,&
  !                                        &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
  !           ENDIF                                       
  !           GEOMETRIC_PARAMETERS(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%&
  !             &PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(GROUP_DERIVATIVES(dev_idx),local_number,0))=NODAL_BUFFER(dev_idx)
  !                          
  !        ENDDO
  !     ENDDO !comp_idx                   
  !     
  !  ENDDO!nn    
  !  
  !  !close a file    
  !  CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)
  !
  !  !release the temporary memory
  !  IF(ALLOCATED(NODAL_BUFFER)) THEN
  !     DEALLOCATE(NODAL_BUFFER)
  !  ENDIF
  !  IF(ALLOCATED(GROUP_DERIVATIVES)) THEN
  !     DEALLOCATE(GROUP_DERIVATIVES)
  !  ENDIF
  !  
  !  CALL EXITS("FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE")
  !  RETURN 1  
  !END SUBROUTINE FIELD_IO_IMPORT_NODES_FROM_LOCAL_FILE  
 
  !
  !================================================================================================================================
  !  

  !>Write all the nodal information from LOCAL_PROCESS_NODAL_INFO_SET to local exnode files       
  SUBROUTINE FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE(LOCAL_PROCESS_NODAL_INFO_SET, NAME, my_computational_node_number, &
    &computational_node_numbers,ERR, ERROR, *)
    !the reason that my_computational_node_number is used in the argument is for future extension
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), INTENT(INOUT):: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN):: computational_node_numbers !<total process number       
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LINE 
    TYPE(VARYING_STRING) :: FILE_NAME !the prefix name of file.  
    TYPE(VARYING_STRING) :: FILE_STATUS !the status for opening file.        
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
    INTEGER(INTG) :: FILE_ID, domain_idx, domain_no, local_number, global_number, MY_DOMAIN_INDEX
    INTEGER(INTG), ALLOCATABLE :: GROUP_DERIVATIVES(:)    
    INTEGER(INTG) :: nn, comp_idx,  dev_idx, NUM_OF_NODAL_DEV, MAX_NUM_OF_NODAL_DERIVATIVES
    REAL(DP), ALLOCATABLE :: NODAL_BUFFER(:)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)

    CALL ENTERS("FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE",ERR,ERROR,*999)    
    
    !get my own computianal node number--be careful the rank of process in the MPI pool 
    !is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
    !will be a secured way to get the number         
    FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(my_computational_node_number,"*",ERR,ERROR))//".exnode"        
    FILE_ID=1029+my_computational_node_number
    MAX_NUM_OF_NODAL_DERIVATIVES=0
    
    IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
       CALL FLAG_ERROR("the nodal information set in input is invalid",ERR,ERROR,*999)            
    ENDIF

    IF(.NOT.ALLOCATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("the nodal information set is not associated with any numbering list",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0) THEN 
       CALL FLAG_ERROR("the nodal information set does not contain any nodes",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(1)%SAME_HEADER==.TRUE.) THEN
       CALL FLAG_ERROR("the first header flag of nodal information set should be false",ERR,ERROR,*999)            
    ENDIF
    
    !NULLIFY(NODAL_BUFFER)    
    !NULLIFY(GROUP_DERIVATIVES)
    
    !open a file 
    FILE_STATUS="REPLACE"
    CALL FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*999)
    
    !write out the group name    
    LINE=FIELD_IO_FILEDS_GROUP_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS, ERR,ERROR)    
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get group namein IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF               
    CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    !write out the number of files    
    LINE=FIELD_IO_MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)      
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get multiple file information in IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF             
    !CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       

    DO nn=1, LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
      
       !tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS         
       global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)
       
       !check whether need to write out the nodal information header  
       IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.) THEN
        !write out the nodal header
           
          CALL FIELD_IO_EXPORT_NODAL_GROUP_HEADER_FORTRAN(LOCAL_PROCESS_NODAL_INFO_SET, nn, MAX_NUM_OF_NODAL_DERIVATIVES, &
          &my_computational_node_number, FILE_ID, ERR,ERROR,*999) 
           !value_idx=value_idx-1 !the len of NODAL_BUFFER                      
          !checking: whether need to allocate temporary memory for Io writing
          IF(ALLOCATED(NODAL_BUFFER)) THEN
             IF(SIZE(NODAL_BUFFER)<MAX_NUM_OF_NODAL_DERIVATIVES) THEN 
                DEALLOCATE(NODAL_BUFFER)
                DEALLOCATE(GROUP_DERIVATIVES)
                ALLOCATE(NODAL_BUFFER(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty nodal buffer in IO writing",ERR,ERROR,*999)
                ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
             ENDIF
          ELSE
             ALLOCATE(NODAL_BUFFER(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty nodal buffer in IO writing",ERR,ERROR,*999)
             ALLOCATE(GROUP_DERIVATIVES(MAX_NUM_OF_NODAL_DERIVATIVES),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty derivative buffer in IO writing",ERR,ERROR,*999)           
          ENDIF
       ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.
                                   
       !write out the components' values of this node in this domain
       DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS   

         !finding the local numbering through the global to local mapping   
          DOMAIN_MAPPING_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
          !get the domain index for this variable component according to my own computional node number
          DO domain_idx=1,DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
             domain_no=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
             IF(domain_no==my_computational_node_number) THEN
                MY_DOMAIN_INDEX=domain_idx
                EXIT !out of loop--domain_idx
             ENDIF
          ENDDO !domain_idX
          local_number=DOMAIN_MAPPING_NODES%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
          !use local domain information find the out the maximum number of derivatives
          DOMAIN_NODES=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
          
          !write out the user numbering of node if comp_idx ==1
          IF(comp_idx==1) THEN 
             LINE=" Node:            "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NODES%NODES(local_number)%USER_NUMBER,"*",ERR,ERROR))
             CALL FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
          ENDIF
          
          !get the nodal partial derivatives
          NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
          GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
          !sort  the partial derivatives
          CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
          DO dev_idx=1, NUM_OF_NODAL_DEV
             NULLIFY(GEOMETRIC_PARAMETERS)
             IF(comp_idx==1) THEN
                CALL FIELD_PARAMETER_SET_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%FIELD,&
                                          &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
             ELSE IF (ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%DOMAIN, &
                &TARGET=LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS (comp_idx-1)%PTR%DOMAIN)) THEN
                CALL FIELD_PARAMETER_SET_GET(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%FIELD,&
                                          &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
             ENDIF                                       
             NODAL_BUFFER(dev_idx)=GEOMETRIC_PARAMETERS(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(comp_idx)%PTR%&
               &PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(GROUP_DERIVATIVES(dev_idx),local_number,0))
             CALL FIELD_IO_FORTRAN_FILE_WRITE_DP(FILE_ID, NODAL_BUFFER, NUM_OF_NODAL_DEV, ERR,ERROR,*999)                
          ENDDO
       ENDDO !comp_idx                   
       
    ENDDO!nn    
    
    !close a file    
    CALL FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)

    !release the temporary memory
    IF(ALLOCATED(NODAL_BUFFER)) DEALLOCATE(NODAL_BUFFER) 
    IF(ALLOCATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)    
    IF(ASSOCIATED(GEOMETRIC_PARAMETERS)) NULLIFY(GEOMETRIC_PARAMETERS)
    
    CALL EXITS("FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE")
    RETURN
999 CALL ERRORS("FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE",ERR,ERROR)
    CALL EXITS("FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE")
    RETURN 1  
  END SUBROUTINE FIELD_IO_EXPORT_NODES_INTO_LOCAL_FILE
    
  !
  !================================================================================================================================
  !  

  !>Get the region label  
  FUNCTION FIELD_IO_FILEDS_GROUP_INFO_GET(FIELDS, ERR, ERROR)
    !Argument variables   
    TYPE(FIELDS_TYPE), INTENT(IN):: FIELDS !<FILEDS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: FIELD_IO_FILEDS_GROUP_INFO_GET !<On exit group name
   
    CALL ENTERS("FIELD_IO_FILEDS_GROUP_INFO_GET",ERR,ERROR,*999)    
   
    IF(ASSOCIATED(FIELDS%REGION)) THEN
       FIELD_IO_FILEDS_GROUP_INFO_GET=" Group name: "//TRIM(FIELDS%REGION%LABEL)
    ELSE
       CALL FLAG_ERROR("fields are not associated with a region, so can not find a group name",ERR,ERROR,*999)     
    ENDIF                
     
    CALL EXITS("FIELD_IO_FILEDS_GROUP_INFO_GET")
    RETURN
999 CALL ERRORS("FIELD_IO_FILEDS_GROUP_INFO_GET",ERR,ERROR)
    CALL EXITS("FIELD_IO_FILEDS_GROUP_INFO_GET")
  END FUNCTION FIELD_IO_FILEDS_GROUP_INFO_GET

  !
  !================================================================================================================================
  !  

  !>Get the number of files      
  FUNCTION FIELD_IO_MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)
    !Argument variables   
    INTEGER(INTG), INTENT(IN) :: computational_node_numbers  
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   !Temporary variables
    TYPE(VARYING_STRING) :: FIELD_IO_MULTI_FILES_INFO_GET !<ON exit multi fiels information
   
    CALL ENTERS("FIELD_IO_MULTI_FILES_INFO_GET",ERR,ERROR,*999)    
   
    FIELD_IO_MULTI_FILES_INFO_GET="files: "//TRIM(NUMBER_TO_VSTRING(computational_node_numbers,"*",ERR,ERROR))
     
    CALL EXITS("FIELD_IO_MULTI_FILES_INFO_GET")
    RETURN
999 CALL ERRORS("FIELD_IO_MULTI_FILES_INFO_GET",ERR,ERROR)
    CALL EXITS("FIELD_IO_MULTI_FILES_INFO_GET")
  END FUNCTION FIELD_IO_MULTI_FILES_INFO_GET  

  !!
  !!================================================================================================================================
  !!  
  !
  !!>Write a string     
  !SUBROUTINE FIELD_IO_MPIIO_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_OF_DATA, ERR,ERROR,*)
  ! 
  !  !Argument variables   
  !  TYPE(VARYING_STRING), INTENT(IN) :: STRING_DATA !<the string data.
  !  INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
  !  INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string    
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(VARYING_STRING) :: TMP_STRING !<temporary string
  !  INTEGER(INTG) :: TMP_LEN, MPI_STATUS    
  !  
  !  CALL ENTERS("FIELD_IO_MPIIO_FILE_WRITE_STRING",ERR,ERROR,*999)
  !      
  !  IF(LEN_OF_DATA==0) THEN
  !     CALL FLAG_ERROR("leng of string is zero",ERR,ERROR,*999) 
  !  ENDIF
  ! 
  !  TMP_STRING=TRIM(STRING_DATA)//"\n"
  !  TMP_LEN=LEN_TRIM(TMP_STRING)
  !  
  !  CALL MPI_FILE_WRITE(FILE_ID, CHAR(TMP_LEN), TMP_LEN, MPI_CHARACTER, MPI_STATUS, ERR)
  !      
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_WRITE_STRING")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_MPIIO_FILE_WRITE_STRING",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_WRITE_STRING")
  !  RETURN 1
  !END SUBROUTINE FIELD_IO_MPIIO_FILE_WRITE_STRING

  !
  !================================================================================================================================
  !  

  !>Read a string using FORTRAN IO     
  SUBROUTINE FIELD_IO_FORTRAN_FILE_READ_STRING(FILE_ID, STRING_DATA, FILE_END, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(INOUT) :: STRING_DATA !<the string data.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    LOGICAL, INTENT(INOUT) :: FILE_END !< file end
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    CHARACTER (LEN=MAXSTRLEN) :: TEMP_STR
    INTEGER(INTG) :: LEN_OF_DATA, IOS 
    
    STRING_DATA=REMOVE(STRING_DATA, 1, LEN(STRING_DATA))
    
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_READ_STRING", ERR, ERROR, *999)
            
    READ(FILE_ID, "(A)", IOSTAT=IOS), TEMP_STR
    
    IF(IOS>=0)  THEN
       FILE_END=.FALSE.
    ELSE
       FILE_END=.TRUE.
    ENDIF
    
    STRING_DATA=TRIM(TEMP_STR)
    LEN_OF_DATA=LEN(STRING_DATA)
    
    IF(LEN_OF_DATA==0) THEN
       CALL FLAG_ERROR("leng of string is zero",ERR,ERROR,*999) 
    ENDIF    
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_READ_STRING")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_READ_STRING",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_READ_STRING")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_READ_STRING

  !
  !================================================================================================================================
  !  

  !>Write a string using FORTRAN IO     
  SUBROUTINE FIELD_IO_FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(IN) :: STRING_DATA !<the string data.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_WRITE_STRING",ERR,ERROR,*999)
        
    IF(LEN_OF_DATA==0) THEN
       CALL FLAG_ERROR("leng of string is zero",ERR,ERROR,*999) 
    ENDIF
   
    WRITE(FILE_ID, "(A)"), CHAR(STRING_DATA) 
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_WRITE_STRING")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_WRITE_STRING",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_WRITE_STRING")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_WRITE_STRING

  !!
  !!================================================================================================================================
  !!  
  !
  !!>Write a real data using MPIIO     
  !SUBROUTINE FIELD_IO_MPIIO_FILE_WRITE_DP(FILE_ID, REAL_DATA, LEN_OF_DATA, ERR,ERROR,*)
  ! 
  !  !Argument variables   
  !  REAL(DP), INTENT(IN) :: REAL_DATA(:) !<the name of file.
  !  INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID    
  !  INTEGER(INTG), INTENT(INOUT) :: LEN_OF_DATA !<length of string    
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(VARYING_STRING) :: STRING_DATA !<the name of file.
  !  INTEGER(INTG) :: idx, MPI_STATUS
  !  
  !  CALL ENTERS("FIELD_IO_MPIIO_FILE_WRITE_DP",ERR,ERROR,*999)    
  ! 
  !  !DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"ES)"
  !  !WRITE(STRING_DATA, CHAR(DP_FMT)), REAL_DATA(1:LEN_OF_DATA)
  !  
  !  DO idx=1,LEN_OF_DATA
  !     STRING_DATA=TRIM(STRING_DATA)//TRIM(NUMBER_TO_VSTRING(REAL_DATA(idx),"*",ERR,ERROR))
  !  ENDDO
  !      
  !  STRING_DATA=TRIM(STRING_DATA)//"\n"
  !  LEN_OF_DATA=LEN_TRIM(STRING_DATA)
  !
  !  CALL MPI_FILE_WRITE(FILE_ID, CHAR(STRING_DATA), LEN_OF_DATA, MPI_CHARACTER, MPI_STATUS, ERR)
  !  
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_WRITE_DP")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_MPIIO_FILE_WRITE_DP",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_WRITE_DP")
  !  RETURN 1
  !END SUBROUTINE FIELD_IO_MPIIO_FILE_WRITE_DP

  !
  !================================================================================================================================
  !  

  !>Read a real data using FORTRAN IO     
  SUBROUTINE FIELD_IO_FORTRAN_FILE_READ_DP(FILE_ID, REAL_DATA, LEN_OF_DATA, FILE_END, ERR,ERROR,*)
  
    !Argument variables   
    REAL(DP), INTENT(OUT) :: REAL_DATA(:) !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    LOGICAL, INTENT(INOUT) :: FILE_END ! file end
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: DP_FMT !<the name of file. 
    INTEGER(INTG) :: IOS   
    
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_READ_DP",ERR,ERROR,*999)    
   
    DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"ES)"
    READ(FILE_ID, CHAR(DP_FMT), IOSTAT=IOS), REAL_DATA(1:LEN_OF_DATA)

    IF(IOS>=0)  THEN
       FILE_END=.FALSE.
    ELSE
       FILE_END=.TRUE.
    ENDIF 
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_READ_DP")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_READ_DP",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_READ_DP")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_READ_DP

  !
  !================================================================================================================================
  !  

  !>Write a real data using FORTRAN IO     
  SUBROUTINE FIELD_IO_FORTRAN_FILE_WRITE_DP(FILE_ID, REAL_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    REAL(DP), INTENT(IN) :: REAL_DATA(:) !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: DP_FMT !<the name of file.    
    
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_WRITE_DP",ERR,ERROR,*999)    
   
    DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"ES)"
    WRITE(FILE_ID, CHAR(DP_FMT)), REAL_DATA(1:LEN_OF_DATA)
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_WRITE_DP")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_WRITE_DP",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_WRITE_DP")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_WRITE_DP

  !!
  !!================================================================================================================================
  !!  
  !
  !!>Write a real data using MPIIO     
  !SUBROUTINE FIELD_IO_MPIIO_FILE_WRITE_INTG(FILE_ID, INTG_DATA, LEN_OF_DATA, ERR,ERROR,*)
  ! 
  !  !Argument variables   
  !  INTEGER(INTG), INTENT(IN) :: INTG_DATA(:) !<the name of file.
  !  INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID    
  !  INTEGER(INTG), INTENT(INOUT) :: LEN_OF_DATA !<length of string    
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  TYPE(VARYING_STRING) :: STRING_DATA !<the name of file.
  !  INTEGER(INTG) :: idx, MPI_STATUS
  !  
  !  CALL ENTERS("FIELD_IO_MPIIO_FILE_WRITE_INTG",ERR,ERROR,*999)    
  ! 
  !  !DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"ES)"
  !  !STRING_DATA=TRIM(NUMBER_TO_VSTRING(INTG_DATA(1:LEN_OF_DATA),"*",ERR,ERROR))
  !  DO idx=1,LEN_OF_DATA
  !     STRING_DATA=TRIM(STRING_DATA)//TRIM(NUMBER_TO_VSTRING(INTG_DATA(idx),"*",ERR,ERROR))
  !  ENDDO
  !  
  !  STRING_DATA=TRIM(STRING_DATA)//"\n"
  !  LEN_OF_DATA=LEN_TRIM(STRING_DATA)
  !
  !  CALL MPI_FILE_WRITE(FILE_ID, CHAR(STRING_DATA), LEN_OF_DATA, MPI_CHARACTER, MPI_STATUS, ERR)
  !  
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_WRITE_INTG")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_MPIIO_FILE_WRITE_INTG",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_WRITE_INTG")
  !  RETURN 1
  !END SUBROUTINE FIELD_IO_MPIIO_FILE_WRITE_INTG

  !
  !================================================================================================================================
  !  

  !>Read a integer data      
  SUBROUTINE FIELD_IO_FORTRAN_FILE_READ_INTG(FILE_ID, INTG_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    INTEGER(INTG), INTENT(OUT) :: INTG_DATA(:) !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: DP_FMT !<the name of file.
       
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_READ_INTG",ERR,ERROR,*999)    
   
    DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"I)"
    READ(FILE_ID, CHAR(DP_FMT)), INTG_DATA(1:LEN_OF_DATA)
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_READ_INTG")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_READ_INTG",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_READ_INTG")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_READ_INTG

  !
  !================================================================================================================================
  !  

  !>Write a integer data      
  SUBROUTINE FIELD_IO_FORTRAN_FILE_WRITE_INTG(FILE_ID, INTG_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    INTEGER(INTG), INTENT(IN) :: INTG_DATA(:) !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: DP_FMT !<the name of file.
       
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_WRITE_INTG",ERR,ERROR,*999)    
   
    DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"I)"
    WRITE(FILE_ID, CHAR(DP_FMT)), INTG_DATA(1:LEN_OF_DATA)
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_WRITE_INTG")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_WRITE_INTG",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_WRITE_INTG")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_WRITE_INTG
 
  !
  !================================================================================================================================
  !  

  !>Open a file using Fortran 
  SUBROUTINE FIELD_IO_FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, FILE_STATUS, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<the name of file.
    TYPE(VARYING_STRING), INTENT(IN) :: FILE_STATUS !<status for opening a file
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
        
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_OPEN",ERR,ERROR,*999)       

    OPEN(UNIT=FILE_ID, FILE=CHAR(FILE_NAME), STATUS=CHAR(FILE_STATUS), FORM="FORMATTED", ERR=999)   
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_OPEN")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_OPEN",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_OPEN")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_OPEN
  
  !!
  !!================================================================================================================================
  !!  
  !
  !!>Open a file using MPIIO 
  !SUBROUTINE FIELD_IO_MPIIO_FILE_OPEN(FILE_ID, FILE_NAME, my_computational_node_number, bufsize, ERR,ERROR,*)
  ! 
  !  !Argument variables   
  !  TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<the name of file.
  !  INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
  !  INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<my computational node number
  !  INTEGER(INTG), INTENT(IN) :: bufsize !<size of buffer 
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  INTEGER(KIND=MPI_OFFSET_KIND) :: DISP
  !      
  !  CALL ENTERS("FIELD_IO_MPIIO_FILE_OPEN",ERR,ERROR,*999)       
  !
  !  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, CHAR(FILE_NAME), MPI_MODE_RDONLY, MPI_INFO_NULL, FILE_ID,  ERR)
  !  DISP=my_computational_node_number*bufsize
  !  CALL MPI_FILE_SEEK(FILE_ID, DISP, MPI_SEEK_SET, ERR) 
  !  
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_OPEN")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_MPIIO_FILE_OPEN",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_OPEN")
  !  RETURN 1
  !END SUBROUTINE FIELD_IO_MPIIO_FILE_OPEN

  !
  !================================================================================================================================
  !  

  !>Close a file using FORTRAN IO to rewind
  SUBROUTINE FIELD_IO_FORTRAN_FILE_REWIND(FILE_ID, ERR,ERROR,*)
  
    !Argument variables   
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_REWIND",ERR,ERROR,*999)    
   
    REWIND(UNIT=FILE_ID, ERR=999)   
    
    CALL EXITS("FIELD_IO_FORTRAN_FILE_REWIND")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_REWIND",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_REWIND")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_REWIND

  !
  !================================================================================================================================
  !
  
  !>Close a file using Fortran
  SUBROUTINE FIELD_IO_FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR, *)
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    CALL ENTERS("FIELD_IO_FORTRAN_FILE_CLOSE",ERR,ERROR,*999)
  
    CLOSE(UNIT=FILE_ID, ERR=999)
  
    CALL EXITS("FIELD_IO_FORTRAN_FILE_CLOSE")
    RETURN
999 CALL ERRORS("FIELD_IO_FORTRAN_FILE_CLOSE",ERR,ERROR)
    CALL EXITS("FIELD_IO_FORTRAN_FILE_CLOSE")
    RETURN 1
  END SUBROUTINE FIELD_IO_FORTRAN_FILE_CLOSE
 
  SUBROUTINE STRING_TO_MUTI_INTEGERS_VS(STRING, NUMBER_OF_INTEGERS, INTG_DATA, ERR, ERROR, *)
  
    !#### Function: STRING_TO_INTEGER_VS
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Converts a varying string representation of a number to an integer.
    !###  Parent-function: STRING_TO_INTEGER
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(INOUT) :: INTG_DATA(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_INTEGERS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING, LOCAL_STRING1
    INTEGER(INTG) :: idx, pos
    
    CALL ENTERS("STRING_TO_MUTI_INTEGERS_VS",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING

    LOCAL_STRING=STRING

    DO idx=1,NUMBER_OF_INTEGERS-1
       LOCAL_STRING=ADJUSTL(LOCAL_STRING)
       LOCAL_STRING=TRIM(LOCAL_STRING)
       pos=INDEX(LOCAL_STRING, " ")
       LOCAL_STRING1=EXTRACT(LOCAL_STRING, 1, pos-1)
       INTG_DATA(idx)=STRING_TO_INTEGER(LOCAL_STRING1, ERR, ERROR)
       LOCAL_STRING=REMOVE(LOCAL_STRING,1,pos)
    ENDDO
    LOCAL_STRING=ADJUSTL(LOCAL_STRING)
    LOCAL_STRING=TRIM(LOCAL_STRING)
    INTG_DATA(idx)=STRING_TO_INTEGER(LOCAL_STRING, ERR, ERROR)
        
    CALL EXITS("STRING_TO_MUTI_INTEGERS_VS")
    RETURN
999 CALL ERRORS("STRING_TO_MUTI_INTEGERS_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_MUTI_INTEGERS_VS")
    RETURN 1
  END SUBROUTINE STRING_TO_MUTI_INTEGERS_VS

  SUBROUTINE STRING_TO_MUTI_REALS_VS(STRING, NUMBER_OF_REALS, REAL_DATA, ERR, ERROR, *)
  
    !#### Function: STRING_TO_INTEGER_VS
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Converts a varying string representation of a number to an integer.
    !###  Parent-function: STRING_TO_INTEGER
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    REAL(DP), INTENT(INOUT) :: REAL_DATA(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_REALS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING, LOCAL_STRING1
    INTEGER(INTG) :: idx, pos
    
    CALL ENTERS("STRING_TO_MUTI_REALS_VS",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING

    LOCAL_STRING=STRING

    DO idx=1,NUMBER_OF_REALS-1
       LOCAL_STRING=ADJUSTL(LOCAL_STRING)
       LOCAL_STRING=TRIM(LOCAL_STRING)
       pos=INDEX(LOCAL_STRING, " ")
       LOCAL_STRING1=EXTRACT(LOCAL_STRING, 1, pos-1)
       REAL_DATA(idx)=STRING_TO_DOUBLE(LOCAL_STRING1, ERR, ERROR)
       LOCAL_STRING=REMOVE(LOCAL_STRING,1,pos)
    ENDDO
    LOCAL_STRING=ADJUSTL(LOCAL_STRING)
    LOCAL_STRING=TRIM(LOCAL_STRING)
    REAL_DATA(idx)=STRING_TO_DOUBLE(LOCAL_STRING, ERR, ERROR)

    CALL EXITS("STRING_TO_MUTI_REALS_VS")
    RETURN
999 CALL ERRORS("STRING_TO_MUTI_REALS_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_MUTI_REALS_VS")
    RETURN 1
  END SUBROUTINE STRING_TO_MUTI_REALS_VS

 
  !!
  !!================================================================================================================================
  !!  
  !
  !!>Close a file using MPIIO 
  !SUBROUTINE FIELD_IO_MPIIO_FILE_CLOSE(FILE_ID, ERR,ERROR,*)
  ! 
  !  !Argument variables   
  !  INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
  !  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
  !  !Local Variables
  !  
  !  CALL ENTERS("FIELD_IO_MPIIO_FILE_CLOSE",ERR,ERROR,*999)    
  ! 
  !  CALL MPI_FILE_CLOSE(FILE_ID, ERR)   
  !  
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_CLOSE")
  !  RETURN
!999 CALL ERRORS("FIELD_IO_MPIIO_FILE_CLOSE",ERR,ERROR)
  !  CALL EXITS("FIELD_IO_MPIIO_FILE_CLOSE")
  !  RETURN 1
  !END SUBROUTINE FIELD_IO_MPIIO_FILE_CLOSE
  
END MODULE FIELD_IO_ROUTINES

   
