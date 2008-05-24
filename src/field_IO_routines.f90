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

  IMPLICIT NONE

  PRIVATE
  
  !Module parameters

  !>Type for lable
  INTEGER(INTG), PARAMETER :: FIELD_LABEL=1  
  INTEGER(INTG), PARAMETER :: VARIABLE_LABEL=2  
  INTEGER(INTG), PARAMETER :: COMPONENT_LABEL=3  
  INTEGER(INTG), PARAMETER :: DERIVATIVE_LABEL=4

  !>Type of scale factor
  INTEGER(INTG), PARAMETER :: SCALE_FACTORS_NUMBER_TYPE=5
  INTEGER(INTG), PARAMETER :: SCALE_FACTORS_PROPERTY_TYPE=6

  !Module types
  
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
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), ALLOCATABLE :: COMPONENTS(:) !<A array of pointers to those components of the node in this local domain       
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
  
  !*******************************************************************************************************************************!
  !for more technonical details about the implementation here, pls go  
  !to read book <<Using MPI2>> written by William Group, Ewing Lusk and
  !Rajeev Thakur. Figures in Page13, Page15 and Page20.
  !    
  !first verion: the data are written out in MPI environment using fortran IO
  !general procedure for exnode: <1> call NODAL_INFO_SET_INITIALISE to initialize a nodal_info_set
  !<2> call IO_INFO_SET_GET to fill the data into nodal_info_set
  !<3> call IO_INFO_SET_SORT to sort the data to group the data with the 
  !same output format and sort the components according to component number 
  !(but this one can be skipped, then every output nodal data will have one header
  !<4> call IO_INFO_SET_EXPORT to export the nodal field variables into exnode files
  !<5> call IO_INFO_SET_FINALIZE to release the IO memory             
  !or just call the simple command EXNODE_EXPORT
  !
  !********************************************************************************************************************************!                 

  !!>export exnode files \see FIELD_IO.
  !INTERFACE FIELDS_NODE_EXPORT
  !  MODULE PROCEDURE EXNODE_INTO_MULTIPLE_FILES 
  !END INTERFACE!FIELDS_NODE_EXPORT

  !!>export exelem files \see FIELD_IO.
  !INTERFACE FIELDS_ELEM_EXPORT
  !  MODULE PROCEDURE EXELEM_INTO_MULTIPLE_FILES 
  !END INTERFACE !FIELDS_ELEM_EXPORT

  !>write data into files \see FIELD_IO
  !INTERFACE CMISS_FILE_WRITE
  !  MODULE PROCEDURE FORTRAN_FILE_WRITE_STRING
  !  MODULE PROCEDURE FORTRAN_FILE_WRITE_DP
  !  MODULE PROCEDURE FORTRAN_FILE_WRITE_INTG
  !END INTERFACE !CMISS_FILE_WRITE     

  !>open a file \see FIELD_IO
  !INTERFACE FORTRAN_FILE_OPEN
  !  !MODULE PROCEDURE MPI_FILE_OPEN
  !  MODULE PROCEDURE FORTRAN_FILE_OPEN
  !END INTERFACE !FORTRAN_FILE_OPEN     

  !!>close a file \see FIELD_IO
  !INTERFACE FORTRAN_FILE_CLOSE
  !  !MODULE PROCEDURE MPI_FILE_CLOSE
  !  MODULE PROCEDURE FORTRAN_FILE_CLOSE
  !END INTERFACE !FORTRAN_FILE_CLOSE      


  !********************************************************************!
  ! Generic-interface: IO_INFO_SET_INITIALISE (public)                 
  ! initialize the nodal information set for IO                        
  !Child-subroutines: NODAL_INFO_SET_INITIALISE(nodal base)            
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !INTERFACE IO_INFO_SET_INITIALISE
  !  MODULE PROCEDURE NODAL_INFO_SET_INITIALISE
  !END INTERFACE !
  !********************************************************************!

  !**************************************************************************************************!
  ! Generic-interface: IO_INFO_SET_GET (public)                                                      
  ! sorting the nodal information set according for IO                                               
  ! Child-subroutines:                                                  
  ! NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(nodal base): collecting the data from each MPI process and filling
  ! into process_nodal_info_set structure    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !INTERFACE IO_INFO_SET_GET
  !  MODULE PROCEDURE NODAL_INFO_SET_ATTACH_LOCAL_PROCESS
  !END INTERFACE !
  !********************************************************************!

  !********************************************************************!
  ! Generic-interface: IO_INFO_SET_SORT (public)                       
  ! sorting the nodal information set for IO                 
  ! Child-subroutines:                                                  
  ! NODAL_INFO_SET_SORT(nodal base): sort according to simality of data,
  ! and also arrange the components according the components number.    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !INTERFACE IO_INFO_SET_SORT
  !  MODULE PROCEDURE NODAL_INFO_SET_SORT
  !END INTERFACE !
  !********************************************************************!

  !***************************************************************************************!
  ! Generic-interface: IO_INFO_SET_EXPORT (public)                       
  ! exporting the fields into disk/files                 
  ! Child-subroutines:                                                  
  ! EXPORT_FIELDS_IN_FORTRAN(nodal base): export the fields into exnode files using Fortran IO    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !INTERFACE IO_INFO_SET_EXPORT
  !  MODULE PROCEDURE EXPORT_FIELDS_IN_NODES_FORTRAN
  !END INTERFACE !
  !***************************************************************************************!


  !********************************************************************!
  ! Generic-interface: IO_INFO_SET_FINALIZE (public)                 
  ! finalize the nodal information set for IO                          
  !Child-subroutines: NODAL_INFO_SET_FINALIZE(nodal base)              
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Interfaces for export node
  !INTERFACE IO_INFO_SET_FINALIZE
  !  MODULE PROCEDURE NODAL_INFO_SET_FINALIZE
  !END INTERFACE !
  !********************************************************************!

  !********************************************************************!
  ! Generic-interface: LABEL_INFO_GET  (private)                       !
  ! get the label for different types                                  !
  !Child-subroutines: LABEL_FIELD_INFO_GET, LABEL_DERIVATIVE_INFO_GET  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Interfaces for export node
  !INTERFACE LABEL_INFO_GET
  !  MODULE PROCEDURE LABEL_FIELD_INFO_GET
  !  MODULE PROCEDURE LABEL_DERIVATIVE_INFO_GET
  !END INTERFACE !
  !********************************************************************!

  PUBLIC :: FIELD_IO_NODES_EXPORT, FIELD_IO_ELEMENTS_EXPORT
 

CONTAINS  
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
    INTEGER(INTG):: my_computational_node_number !<local process number    
    INTEGER(INTG):: computational_node_numbers   !<total process number      
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), ALLOCATABLE:: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<elemental information in this process

    CALL ENTERS("EXELEM_INTO_MULTIPLE_FILES", ERR,ERROR,*999)    
       
    !Get the number of computational nodes
    computational_node_numbers=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !Get my computational node number
    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999

    IF(STRCMP(METHOD, "FORTRAIN") THEN
       ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET),STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate elemental buffer in IO writing",ERR,ERROR,*999)

       CALL ELEMENTAL_INFO_SET_INITIALISE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
       CALL ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*999)
       CALL ELEMENTAL_INFO_SET_SORT(LOCAL_PROCESS_ELEMENTAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)    
       CALL EXPORT_ELEMENTS_INTO_LOCAL_FILE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FILE_NAME, my_computational_node_number, computational_node_numbers, &
            &ERR, ERROR, *999)
       CALL ELEMENTAL_INFO_SET_FINALIZE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*999)
    ELSE IF(STRCMP(METHOD, "MPIIO") THEN
       CALL FLAG_ERROR("MPIIO is not implemented yet",ERR,ERROR,*999)
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
  FUNCTION BASIS_LHTP_FAMILY_LABEL(BASIS, num_scl, num_node, LABEL_TYPE, ERR, ERROR)
    !Argument variables   
    TYPE(BASIS_TYPE), INTENT(IN) :: BASIS !<The error string
    INTEGER(INTG), INTENT(INOUT) :: num_scl
    INTEGER(INTG), INTENT(INOUT) :: num_node
    INTEGER(INTG), INTENT(IN) ::LABEL_TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni 
    TYPE(VARYING_STRING) ::BASIS_LHTP_FAMILY_LABEL
   
    CALL ENTERS("BASIS_LHTP_FAMILY_LABEL",ERR,ERROR,*999)    
   
    IF(BASIS%NUMBER_OF_XI==0) CALL FLAG_ERROR("number of xi in the basis is zero",ERR,ERROR,*999)           
    
    num_scl=1;
    num_node=1
    DO ni=1,BASIS%NUMBER_OF_XI-1
       SELECT CASE(BASIS%INTERPOLATION_XI(ni))
         CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"l.Lagrange*"
             num_scl=num_scl*2  
             num_node=num_node*2
          CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"q.Lagrange*"
             num_scl=num_scl*3
             num_node=num_node*3
          CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"c.Lagrange*"
             num_scl=num_scl*4
             num_node=num_node*4
          CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"c.Hermite*"
             num_scl=num_scl*2*2
             num_node=num_node*2
          CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"q1.Hermite*"
             num_scl=num_scl*2*2
             num_node=num_node*2            
          CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"q2.Hermite*"
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
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"l.Lagrange"
             num_scl=num_scl*2  
             num_node=num_node*2
          CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"q.Lagrange"
             num_scl=num_scl*3
             num_node=num_node*3
          CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"c.Lagrange"
             num_scl=num_scl*4
             num_node=num_node*4
          CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"c.Hermite"
             num_scl=num_scl*2*2
             num_node=num_node*2
          CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"q1.Hermite"
             num_scl=num_scl*2*2
             num_node=num_node*2            
          CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
             BASIS_LHTP_FAMILY_LABEL=BASIS_LHTP_FAMILY_LABEL//"q2.Hermite"
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

    !BASIS_LHTP_FAMILY_LABEL=" "!BASIS_LHTP_FAMILY_LABEL(1:(LEN_TRIM(BASIS_LHTP_FAMILY_LABEL)-1))
    SELECT CASE(LABEL_TYPE)
       CASE (SCALE_FACTORS_NUMBER_TYPE)
          BASIS_LHTP_FAMILY_LABEL=TRIM(BASIS_LHTP_FAMILY_LABEL)//", #Scale factor="//TRIM(NUMBER_TO_VSTRING(num_scl,"*",ERR,ERROR))
       CASE (SCALE_FACTORS_PROPERTY_TYPE)
          BASIS_LHTP_FAMILY_LABEL=TRIM(BASIS_LHTP_FAMILY_LABEL)//", no modify, standard node based"
       CASE DEFAULT
          CALL FLAG_ERROR("Invalid interpolation type",ERR,ERROR,*999)
    END SELECT   
    
    !MAX_SCALE_FACTORS=MAX(MAX_SCALE_FACTORS,num_scl)
                     
    CALL EXITS("BASIS_LHTP_FAMILY_LABEL")
    RETURN
999 CALL ERRORS("BASIS_LHTP_FAMILY_LABEL",ERR,ERROR)
    CALL EXITS("BASIS_LHTP_FAMILY_LABEL")
  END FUNCTION BASIS_LHTP_FAMILY_LABEL    

  !
  !================================================================================================================================
  !  

  !>Write the header of a group elements 
  SUBROUTINE EXPORT_ELEMENTAL_GROUP_HEADER(LOCAL_PROCESS_ELEMENTAL_INFO_SET, LOCAL_ELEMENTAL_NUMBER,  MAX_NODE_CPMP_INDEX, &
     &   NUM_OF_SCALING_FACTOR_SETS, LIST_COMP_SCALE, my_computational_node_number, FILE_ID, ERR,ERROR, *)
    !Argument variables  
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), POINTER :: LOCAL_PROCESS_ELEMENTAL_INFO_SET  !<LOCAL_PROCESS_NODAL_INFO_SET           
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENTAL_NUMBER !<element number
    INTEGER(INTG), INTENT(INOUT) ::  MAX_NODE_CPMP_INDEX !<MAX_NODE_INDEX
    INTEGER(INTG), INTENT(INOUT) :: NUM_OF_SCALING_FACTOR_SETS !<NUM_OF_SCALING_FACTOR_SETS
    INTEGER(INTG), INTENT(INOUT) :: LIST_COMP_SCALE(:)
    INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN) :: FILE_ID !< FILE ID    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    !INTEGER(INTG), INTENT(INOUT) :: MAX_SCALE_FACTORS !<MAX_SCALE_FACTORS
    INTEGER(INTG) :: NUM_OF_VARIABLES, MAX_NUM_NODES !NUM_OF_FIELDS NUM_OF_NODES
    INTEGER(INTG) :: domain_idx, domain_no, MY_DOMAIN_INDEX, local_number, global_number
    INTEGER(INTG), ALLOCATABLE :: LIST_SCALE_FACTORS(:),  GROUP_LOCAL_NUMBER(:), GROUP_SCALE_FACTORS(:), GROUP_NODE(:), GROUP_VARIABLES(:)
    !INTEGER(INTG), ALLOCATABLE :: GROUP_DERIVATIVES(:)
    INTEGER(INTG) ::nn, np, ny2, mk, nk, num_scl, num_node, comp_idx, comp_idx1, scl_idx, scl_idx1, var_idx!value_idx field_idx global_var_idx   
    LOGICAL :: SWITCH   
    TYPE(FIELD_TYPE), POINTER :: field_ptr
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable_ptr   
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER :: tmp_components(:)        
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_ELEMENTS !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS ! domain nodes
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
    TYPE(BASIS_TYPE), POINTER :: BASIS 
    TYPE(BASIS_TYPE), POINTER :: BASIS1
    TYPE(VARYING_STRING) :: LINE, LABEL
   
    CALL ENTERS("EXPORT_ELEMENTAL_GROUP_HEADER",ERR,ERROR,*999)    
    
    !colllect nodal header information for IO first
    
    !!get the number of this computational node from mpi pool
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
    !IF(ERR/=0) GOTO 999        
    
    !attach the temporary pointer
    tmp_components=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS    
    
    !collect maximum number of nodal derivatives, number of fields and variables 
    NUM_OF_SCALING_FACTOR_SETS=0
    NUM_OF_VARIABLES=0
    MAX_NUM_NODES=0
    MAX_NODE_CPMP_INDEX=0
    global_number=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_ELEMENTAL_NUMBER)
    NULLIFY(field_ptr)     
    NULLIFY(variable_ptr)
    IF(ASSOCIATED(GROUP_LOCAL_NUMBER)) DEALLOCATE(GROUP_LOCAL_NUMBER)
    ALLOCATE(GROUP_LOCAL_NUMBER(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate GROUP_LOCAL_NUMBER in exelem header",ERR,ERROR,*999)
    IF(ASSOCIATED(LIST_SCALE_FACTORS)) DEALLOCATE(LIST_SCALE_FACTORS)
    ALLOCATE(LIST_SCALE_FACTORS(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate LIST_SCALE_FACTORS in exelem header",ERR,ERROR,*999)
     
    DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of fields
       !IF (.NOT.ASSOCIATED(field_ptr, target=tmp_components(comp_idx)%PTR%FIELD)) THEN
       !      NUM_OF_FIELDS=NUM_OF_FIELDS+1
       !    field_ptr=>tmp_components(comp_idx)%PTR%FIELD
       !ENDIF         
       
       !calculate the number of variables
       IF (.NOT.ASSOCIATED(variable_ptr, target=tmp_components(comp_idx)%PTR%FIELD_VARIABLE)) THEN
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          variable_ptr=>tmp_components(comp_idx)%PTR%FIELD_VARIABLE
       ENDIF      
              
       !finding the local numbering through the global to local mapping   
       DOMAIN_MAPPING_ELEMENTS=>tmp_components(comp_idx)%PTR%DOMAIN%MAPPINGS%ELEMENTS       
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
       DOMAIN_ELEMENTS=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%ELEMENTS 
       DOMAIN_NODES=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%COMPONENTS(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES
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
                BASIS1=>tmp_components(LIST_SCALE_FACTORS(scl_idx1))%PTR%DOMAIN%TOPOLOGY%ELEMENTS% &
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
    IF(ASSOCIATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)    
    ALLOCATE(GROUP_VARIABLES(NUM_OF_VARIABLES),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
    
    !!Allocate the momery for group of maximum number of derivatives     
    IF(ASSOCIATED(GROUP_SCALE_FACTORS)) DEALLOCATE(GROUP_SCALE_FACTORS)    
    ALLOCATE(GROUP_SCALE_FACTORS(NUM_OF_SCALING_FACTOR_SETS),STAT=ERR)    
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)
    
    IF(ASSOCIATED(GROUP_NODE)) DEALLOCATE(GROUP_NODE)    
    ALLOCATE(GROUP_NODE(NUM_OF_SCALING_FACTOR_SETS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporaty variable buffer in IO",ERR,ERROR,*999)            
        
    !fill information into the group of fields and variables
    NULLIFY(variable_ptr) 
    GROUP_VARIABLES(:)=0         
    NUM_OF_VARIABLES=0
    DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(LOCAL_ELEMENTAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of variables
       IF (.NOT.ASSOCIATED(variable_ptr, target=tmp_components(comp_idx)%PTR%FIELD_VARIABLE)) THEN
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          variable_ptr=>tmp_components(comp_idx)%PTR%FIELD_VARIABLE
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
       ELSE
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1   
       ENDIF                    
    ENDDO  !comp_idx                
    
    !write out the scale factor set information
    LINE="Scale factor sets ="//TRIM(NUMBER_TO_VSTRING(NUM_OF_SCALING_FACTOR_SETS,"*",ERR,ERROR))
    CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    !MAX_SCALE_FACTORS=0
    GROUP_SCALE_FACTORS(:)=0
    DO scl_idx=1,NUM_OF_SCALING_FACTOR_SETS                       
       BASIS=>tmp_components(LIST_SCALE_FACTORS(scl_idx))%PTR%DOMAIN%TOPOLOGY%ELEMENTS% &
         & ELEMENTS(GROUP_LOCAL_NUMBER(LIST_SCALE_FACTORS(scl_idx)))%BASIS
       num_scl=0
       IF(ASSOCIATED(BASIS)) THEN
          SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
             LABEL=BASIS_LHTP_FAMILY_LABEL(BASIS, num_scl, num_node, SCALE_FACTORS_NUMBER_TYPE, ERR,ERROR)
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
       LINE=LABEL
       CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                                    
    ENDDO!scl_idx

    LINE="#NODE=     "//TRIM(NUMBER_TO_VSTRING(MAX_NUM_NODES,"*",ERR,ERROR))
    CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)      
    LINE="#FIELD="//TRIM(NUMBER_TO_VSTRING(NUM_OF_VARIABLES,"*",ERR,ERROR))
    CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)      
        
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
       !IF (.NOT.ASSOCIATED(field_ptr, target=tmp_components(comp_idx)%PTR%FIELD)) THEN
       !      NUM_OF_FIELDS=NUM_OF_FIELDS+1
       !    field_ptr=>tmp_components(comp_idx)%PTR%FIELD
       !ENDIF         
       
       !calculate the number of variables
       !grouping field variables and components together
       IF(.NOT.ASSOCIATED(variable_ptr,TARGET=tmp_components(comp_idx)%PTR%FIELD_VARIABLE)) THEN !different variables            
          var_idx=var_idx+1
          variable_ptr=>tmp_components(comp_idx)%PTR%FIELD_VARIABLE          
          !write out the field information
          LABEL="  "//TRIM(NUMBER_TO_VSTRING(var_idx,"*",ERR,ERROR))//") "&
          &//LABEL_FIELD_INFO_GET(tmp_components(comp_idx)%PTR, VARIABLE_LABEL,ERR,ERROR)
          IF(ERR/=0) THEN
             CALL FLAG_ERROR("can not get variable label",ERR,ERROR,*999)     
             GOTO 999               
          ENDIF        
          LINE=TRIM(LABEL)//", #Components="//TRIM(NUMBER_TO_VSTRING(GROUP_VARIABLES(var_idx),"*",ERR,ERROR))                  
          CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
       ENDIF

       !write out the component information
       LABEL="   "//LABEL_FIELD_INFO_GET(tmp_components(comp_idx)%PTR, COMPONENT_LABEL,ERR,ERROR)
       IF(ERR/=0) THEN
          CALL FLAG_ERROR("can not get component label",ERR,ERROR,*999)     
          GOTO 999               
       ENDIF        
       LINE=TRIM(LABEL)//"."                          
       BASIS=>tmp_components(comp_idx)%PTR%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(comp_idx)%BASIS
       SELECT CASE(BASIS%TYPE)
         CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            LABEL=BASIS_LHTP_FAMILY_LABEL(BASIS, num_scl, num_node, SCALE_FACTORS_PROPERTY_TYPE, ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get basis type of lagrange_hermite label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF                           
            !CASE(BASIS_SIMPLEX_TYPE)
            !  CALL BASIS_SIMPLEX_FAMILY_CREATE(BASIS,ERR,ERROR,*999)
         CASE DEFAULT
            CALL FLAG_ERROR("Basis type "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))//" is invalid or not implemented",ERR,ERROR,*999)
       END SELECT
       LINE=TRIM(LINE)//"  "//TRIM(LABEL) 
       CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999) 
       
       LINE="     #NODE=     "//TRIM(NUMBER_TO_VSTRING(num_node,"*",ERR,ERROR))
       CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)      

       IF(BASIS%DEGENERATE==.FALSE.) THEN
          IF(LIST_COMP_SCALE(comp_idx)==1) THEN
             scl_idx=0
          ELSE   
             scl_idx= SUM(GROUP_SCALE_FACTORS(1:LIST_COMP_SCALE(comp_idx)-1))
          ENDIF   
          BASIS=>DOMAIN_ELEMENTS%ELEMENTS(LOCAL_ELEMENTAL_NUMBER)%BASIS    
          DO nn=1,BASIS%NUMBER_OF_NODES
             LINE="     "//TRIM(NUMBER_TO_VSTRING(nn,"*",ERR,ERROR))//". #Values=     "//&
                 &TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_DERIVATIVES(nn),"*",ERR,ERROR))
             CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                     
             np=DOMAIN_ELEMENTS%ELEMENTS(LOCAL_ELEMENTAL_NUMBER)%ELEMENT_NODES(nn)               
             LINE="      Value indices:     "             
             DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                nk=DOMAIN_ELEMENTS%ELEMENTS(LOCAL_ELEMENTAL_NUMBER)%ELEMENT_DERIVATIVES(mk,nn)
                !ny2=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk)
                !GROUP_DERIVATIVES(mk)=ny2
                LINE=LINE//TRIM(NUMBER_TO_VSTRING(nk, "*",ERR,ERROR))
             ENDDO !mk
             CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                 
  
             LINE="     Scale factor indices:     "             
             DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                scl_idx=scl_idx+1
                LINE=LINE//TRIM(NUMBER_TO_VSTRING(scl_idx, "*",ERR,ERROR))
             ENDDO !mk
               CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                                         
          ENDDO !nn
       ELSE
          CALL FLAG_ERROR("exporting degenerated nodes has not been implemented",ERR,ERROR,*999) 
       ENDIF  
    ENDDO !comp_idx         
    
    
    !release temporary memory
    IF(ASSOCIATED(LIST_SCALE_FACTORS)) DEALLOCATE(LIST_SCALE_FACTORS)
    IF(ASSOCIATED(GROUP_LOCAL_NUMBER)) DEALLOCATE(GROUP_LOCAL_NUMBER)
    IF(ASSOCIATED(GROUP_SCALE_FACTORS)) DEALLOCATE(GROUP_SCALE_FACTORS)    
    IF(ASSOCIATED(GROUP_NODE)) DEALLOCATE(GROUP_NODE)    
    IF(ASSOCIATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)
    !IF(ASSOCIATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)
        
    CALL EXITS("EXPORT_ELEMENTAL_GROUP_HEADER")
    RETURN
999 CALL ERRORS("EXPORT_ELEMENTAL_GROUP_HEADER",ERR,ERROR)
    CALL EXITS("EXPORT_ELEMENTAL_GROUP_HEADER")
    RETURN 1       
  END SUBROUTINE EXPORT_ELEMENTAL_GROUP_HEADER  
  
  !
  !================================================================================================================================
  !  

  !>Write all the elemental information from LOCAL_PROCESS_NODAL_INFO_SET to exelem files 
  SUBROUTINE EXPORT_ELEMENTS_INTO_LOCAL_FILE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, NAME, my_computational_node_number, &
  &computational_node_numbers,ERR, ERROR, *)
    !the reason that my_computational_node_number is used in the argument is for future extension
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), POINTER :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<nodal information in this process
    TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN):: computational_node_numbers !<total process number       
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: FILE_NAME, LINE !the prefix name of file.    
    TYPE(BASIS_TYPE), POINTER ::BASIS
    INTEGER(INTG) :: FILE_ID, domain_idx, domain_no, local_number, global_number, MY_DOMAIN_INDEX, MAX_NODE_CPMP_INDEX
    INTEGER(INTG), POINTER :: LIST_COMP_SCALE(:)!LIST_COMP(:) !Components which will be used for export scale factors    
    INTEGER(INTG) :: nk, np, nn, ns, mk, ny2, elem_num, elem_idx, comp_idx,  scal_idx, NUM_OF_SCALING_FACTOR_SETS !dev_idx   
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_ELEMENTS !The domain mapping to calculate elemental mappings
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS ! domain elements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain elements    
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER :: tmp_components(:)            
    REAL(DP), ALLOCATABLE :: SCALE_FACTORS(:), SCALING_BUFFER(:) 

    CALL ENTERS("EXPORT_ELEMENTS_INTO_LOCAL_FILE",ERR,ERROR,*999)    
    
    !get my own computianal node number--be careful the rank of process in the MPI pool 
    !is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
    !will be a secured way to get the number         
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)         
    !IF(ERR/=0) GOTO 999
    FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(my_computational_node_number,"*",ERR,ERROR))//".exelem"        
    FILE_ID=1030
    NUM_OF_SCALING_FACTOR_SETS=0
    
    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)) THEN
       CALL FLAG_ERROR("the elemental information set in input is invalid",ERR,ERROR,*999)            
    ENDIF

    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("the elemental information set is not associated with any numbering list",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS==0) THEN 
       CALL FLAG_ERROR("the elemental information set does not contain any nodes",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(1)%SAME_HEADER==.TRUE.) THEN
       CALL FLAG_ERROR("the first header flag of elemental information set should be false",ERR,ERROR,*999)            
    ENDIF
    
    NULLIFY(SCALE_FACTORS)    
    NULLIFY(SCALING_BUFFER)
    NULLIFY(LIST_COMP_SCALE)
    NULLIFY(tmp_components)
    
    !open a file 
    CALL FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, ERR,ERROR,*999)
    
    !write out the group name    
    LINE=FILEDS_GROUP_INFO_GET(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS, ERR,ERROR)    
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get group name in IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF               
    CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    !write out the number of files    
    LINE=MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)      
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get multiple file information in IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF             
    !CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       

    ns=1!
    DO elem_idx=1, LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
      
      tmp_components=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%COMPONENTS         
       global_number=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(elem_idx)
       

       IF(ASSOCIATED(LIST_COMP_SCALE)) DEALLOCATE(LIST_COMP_SCALE)         
       ALLOCATE(LIST_COMP_SCALE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%NUMBER_OF_COMPONENTS),STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate LIST_COMP_SCALE in exelem io",ERR,ERROR,*999)

       !check whether need to write out the nodal information header  
       IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%SAME_HEADER==.FALSE.) THEN
          !print * ,"elem_idx:", elem_idx, "in", my_computational_node_number
           !write out the nodal header
           CALL EXPORT_ELEMENTAL_GROUP_HEADER(LOCAL_PROCESS_ELEMENTAL_INFO_SET, elem_idx, MAX_NODE_CPMP_INDEX, NUM_OF_SCALING_FACTOR_SETS, &
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
       LINE="   Element:     "//TRIM(NUMBER_TO_VSTRING(tmp_components(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(tmp_components(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%USER_NUMBER,"*",ERR,ERROR))//" 0 0"
       CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
       !node info
       LINE="   Nodes:" !NEED TO LIST ALL THE NODES(SUCH CUBIC CASE)
       BASIS=>tmp_components(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(tmp_components(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%BASIS
       IF(BASIS%DEGENERATE==.FALSE.) THEN
          DO nn=1,BASIS%NUMBER_OF_NODES
             LINE=LINE//" "//TRIM(NUMBER_TO_VSTRING(tmp_components(MAX_NODE_CPMP_INDEX)%PTR%DOMAIN%MESH%TOPOLOGY(tmp_components(MAX_NODE_CPMP_INDEX)%PTR%MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_number)%USER_ELEMENT_NODES(nn),"*",ERR,ERROR))    
          ENDDO !nn
       ELSE
          CALL FLAG_ERROR("exporting degenerated nodes has not been implemented",ERR,ERROR,*999) 
       ENDIF        
       CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)
                         
       !write out scale factors information        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       LINE="   Scale factor:"
       CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       scal_idx=1
       DO comp_idx=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(elem_idx)%NUMBER_OF_COMPONENTS          

         !finding the local numbering through the global to local mapping   
          DOMAIN_MAPPING_ELEMENTS=>tmp_components(comp_idx)%PTR%DOMAIN%MAPPINGS%ELEMENTS       
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
          DOMAIN_ELEMENTS=>tmp_components(comp_idx)%PTR%DOMAIN%TOPOLOGY%ELEMENTS
         DOMAIN_NODES=>tmp_components(comp_idx)%PTR%DOMAIN%TOPOLOGY%NODES      
                    
          !write out the components' values of this node in this domain
          !DO scal_idx=1, NUM_OF_SCALING_FACTOR_SETS   
          IF(LIST_COMP_SCALE(comp_idx)==scal_idx) THEN
             scal_idx=scal_idx+1
             !global_number=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(elem_num)

            !finding the local numbering through the global to local mapping   
             !DOMAIN_MAPPING_ELEMENTS=>tmp_components(LIST_COMP(scal_idx))%PTR%DOMAIN%MAPPINGS%ELEMENTS       
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
             IF(ASSOCIATED(SCALING_BUFFER)) DEALLOCATE(SCALING_BUFFER)            
             ALLOCATE(SCALING_BUFFER(SUM(BASIS%NUMBER_OF_DERIVATIVES(1:BASIS%NUMBER_OF_NODES))),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate scale buffer in IO",ERR,ERROR,*999)
          
             BASIS=>DOMAIN_ELEMENTS%ELEMENTS(local_number)%BASIS   
                           
             !CALL DISTRIBUTED_VECTOR_DATA_GET(tmp_components(comp_idx)%PTR%FIELD%SCALINGS%SCALINGS(tmp_components(comp_idx)%PTR% &
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
             CALL FORTRAN_FILE_WRITE_DP(FILE_ID, SCALING_BUFFER, ns, ERR,ERROR,*999)
             !IF(ASSOCIATED(SCALING_BUFFER)) DEALLOCATE(SCALING_BUFFER)              
          ENDIF !(LIST_COMP_SCALE(comp_idx)==scal_idx)LIST_COMP(:)
       ENDDO !comp_idx=1         
    ENDDO!elem_idx    
    
    !close a file    
    CALL FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)

    !release the temporary memory
    IF(ASSOCIATED(SCALING_BUFFER)) DEALLOCATE(SCALING_BUFFER)
    IF(ASSOCIATED(SCALE_FACTORS)) NULLIFY(SCALE_FACTORS)
    !IF(ASSOCIATED(LIST_COMP)) DEALLOCATE(LIST_COMP)
    IF(ASSOCIATED(LIST_COMP_SCALE)) DEALLOCATE(LIST_COMP_SCALE)
    
    
    CALL EXITS("EXPORT_ELEMENTS_INTO_LOCAL_FILE")
    RETURN
999 CALL ERRORS("EXPORT_ELEMENTS_INTO_LOCAL_FILE",ERR,ERROR)
    CALL EXITS("EXPORT_ELEMENTS_INTO_LOCAL_FILE")
    RETURN 1  
  END SUBROUTINE EXPORT_ELEMENTS_INTO_LOCAL_FILE

  !
  !================================================================================================================================
  !  

  !>Sort the Elemental_info_set according to the type of field variable components 
  SUBROUTINE ELEMENTAL_INFO_SET_SORT(LOCAL_PROCESS_ELEMENTAL_INFO_SET, my_computational_node_number, ERR,ERROR,*)      
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), POINTER :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<elemental information in this process
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: domain_idx,domain_no, MY_DOMAIN_INDEX !temporary variable
    INTEGER(INTG) :: global_number1, local_number1, global_number2, local_number2
    INTEGER(INTG) :: component_idx, nn1, nn2 ! nn, tmp2, tmp1!temporary variable
    LOGICAL :: SWITCH
    !INTEGER(INTG), ALLOCATABLE:: LIST_OF_GLOBAL_NUMBER(:)
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER :: tmp_components(:)        
    !TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: tmp_ptr !temporary variable component    
    !TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !field variable
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_ELEMENTS !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS1, DOMAIN_ELEMENTS2! domain nodes

    !from now on, global numbering are used
    CALL ENTERS("ELEMENTAL_INFO_SET_SORT",ERR,ERROR,*999)    
  
    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("list of global numbering in the input data is invalid",ERR,ERROR,*999)     
    ENDIF
    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)) THEN
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
                  DOMAIN_ELEMENTS1=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%COMPONENTS(component_idx)%PTR%DOMAIN%TOPOLOGY%ELEMENTS
               
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
                  DOMAIN_ELEMENTS2=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS(component_idx)%PTR%DOMAIN%TOPOLOGY%ELEMENTS                   
                   
                  !checking whether they have the same basis
                  IF(DOMAIN_ELEMENTS1%ELEMENTS(local_number1)%BASIS%GLOBAL_NUMBER/=&
                     &DOMAIN_ELEMENTS2%ELEMENTS(local_number2)%BASIS%GLOBAL_NUMBER) THEN                    
                     SWITCH=.FALSE.
                     EXIT     
                  ENDIF   !DOMAIN_ELEMENTS1
               ENDDO !component_idx
            ENDIF !SWITCH==.TRUE.
         ENDIF !LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS==LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn+1)%NUMBER_OF_COMPONENTS
          
          !find two nodes which have the same output, and then they should put together
          IF(SWITCH==.TRUE.) THEN
             !exchange the pointer(to the list the components)
             tmp_components=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%COMPONENTS=>&
             &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%COMPONENTS
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%COMPONENTS=>tmp_components
                          
             !setting the header information
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%SAME_HEADER=.FALSE.
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%SAME_HEADER=.TRUE.
             
             !exchange the number of components
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn2)%NUMBER_OF_COMPONENTS=&
             &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1+1)%NUMBER_OF_COMPONENTS=&
             &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn1)%NUMBER_OF_COMPONENTS

             !exchange the global number
             LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn2)=LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn1+1)
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

    CALL EXITS("ELEMENTAL_INFO_SET_SORT")
    RETURN
999 CALL ERRORS("ELEMENTAL_INFO_SET_SORT",ERR,ERROR)
    CALL EXITS("ELEMENTAL_INFO_SET_SORT")
    RETURN 1  
  END SUBROUTINE ELEMENTAL_INFO_SET_SORT
  
  !
  !================================================================================================================================
  !  

  !>Collect the elemental information from each MPI process     
  SUBROUTINE ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), POINTER:: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT):: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx, var_idx, component_idx, np, nn!temporary variable
    LOGICAL :: SWITCH
    INTEGER(INTG), ALLOCATABLE:: NEW_LIST(:)
    TYPE(FIELD_TYPE), POINTER :: FIELD
    !TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !loca domain
    TYPE(DOMAIN_MAPPING_TYPE), POINTER:: DOMAIN_ELEMENTS_MAPPING !nodes in local mapping--it is different as exnode
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER:: tmp_comp(:) !component of field variable
    TYPE(FIELD_VARIABLE_TYPE), POINTER:: FIELD_VARIABLE !field variable

    CALL ENTERS("ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS",ERR,ERROR,*999)    
   
    !initialize the pointer
    NULLIFY(NEW_LIST)
   
    !attache local process to local nodal information set. In current opencmiss system,
    !each local process owns it local nodal information, so all we need to do is to fill the nodal 
    !information set with nodal information of local process 
    IF((LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS==0).AND.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS) &
    & .AND.(.NOT.ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET))) THEN  
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
                        IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(np)) &
                        &THEN 
                           SWITCH=.TRUE.
                           EXIT
                        ENDIF                    
                     ENDDO
                     !have one more global node
                     !i hate the codes here, but i have to save the memory 
                     IF(SWITCH==.FALSE.) THEN
                        IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS/=0) THEN
                           !crease a new memory space
                           ALLOCATE(NEW_LIST(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS+1),STAT=ERR)
                             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                             !add one more node
                           NEW_LIST(1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS)=&
                           &LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS)              
                           NEW_LIST(LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS+1)=&
                           &DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(np)
                           !release the old memory space
                           DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>NEW_LIST
                           NULLIFY(NEW_LIST)
                           LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS+1
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
         NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS)                                  
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
                       IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)==DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(np)) &
                       &THEN 
                          EXIT
                       ENDIF                    
                    ENDDO
                    !allocate variable component memory
                    IF(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS/=0) THEN
                        tmp_comp=>LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS
                        NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS)
                       
                           ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(&
                           &LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                 
                     LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(&
                     &1:LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS)=tmp_comp(:)
                 
                          LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%&
                          &NUMBER_OF_COMPONENTS+1)%PTR=>FIELD_VARIABLE%COMPONENTS(component_idx)                                       
                     ELSE
                           ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                          LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nn)%COMPONENTS(1)%PTR=>FIELD_VARIABLE%COMPONENTS(component_idx)
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
      NULLIFY(tmp_comp)
     
      !LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>LIST_OF_GLOBAL_NUMBER
      !NULLIFY(LIST_OF_GLOBAL_NUMBER)      
      !release the temporary meomery           
      IF(ASSOCIATED(NEW_LIST)) DEALLOCATE(NEW_LIST)             
    ELSE
      CALL FLAG_ERROR("nodal information set is not initialized properly, and call start method first",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN
999 CALL ERRORS("ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS",ERR,ERROR)
    CALL EXITS("ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN 1  
  END SUBROUTINE ELEMENTAL_INFO_SET_ATTACH_LOCAL_PROCESS      

  !
  !================================================================================================================================
  !  

  !>Finalized the elemental information set     
  SUBROUTINE ELEMENTAL_INFO_SET_FINALIZE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), POINTER :: LOCAL_PROCESS_ELEMENTAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nele, ncomp  !temporary variable

    CALL ENTERS("ELEMENTAL_INFO_SET_FINALIZE",ERR,ERROR,*999)    
   
    IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)) THEN
       NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)
    ENDIF
    DO nele=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
       IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)) THEN
          DO ncomp=1,LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%NUMBER_OF_COMPONENTS
             IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS(ncomp)%PTR)) THEN
                NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS(ncomp)%PTR)
             ENDIF
          ENDDO
          LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%NUMBER_OF_COMPONENTS=0
          LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%SAME_HEADER=.FALSE.
          IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)) THEN
             DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)
          ENDIF   
       ENDIF
    ENDDO
   
    LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=0
    IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    ENDIF
            
    CALL EXITS("ELEMENTAL_INFO_SET_FINALIZE")
    RETURN
999 CALL ERRORS("ELEMENTAL_INFO_SET_FINALIZE",ERR,ERROR)
    CALL EXITS("ELEMENTAL_INFO_SET_FINALIZE")
    RETURN 1  
  END SUBROUTINE ELEMENTAL_INFO_SET_FINALIZE

  !
  !================================================================================================================================
  !  

  !>Initialize the elemental information set         
  SUBROUTINE ELEMENTAL_INFO_SET_INITIALISE(LOCAL_LOCAL_PROCESS_ELEMENTAL_INFO_SET, FIELDS, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_ELEMENTALL_INFO_SET), POINTER :: LOCAL_LOCAL_PROCESS_ELEMENTAL_INFO_SET !<elemental information in this process
    TYPE(FIELDS_TYPE), POINTER ::FIELDS !<the field object
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nfid, nele, ncomp  !temporary variable
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELEMENTAL_INFO_SET_INITIALISE",ERR,ERROR,*999)    

    !validate the input data
    !checking whether the list of fields in the same region
    IF(ASSOCIATED(FIELDS%REGION)) THEN
       DO nfid =2, FIELDS%NUMBER_OF_FIELDS
          IF(FIELDS%FIELDS(nfid-1)%PTR%REGION%USER_NUMBER/=FIELDS%FIELDS(nfid)%PTR%REGION%USER_NUMBER) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(nfid-1,"*",ERR,ERROR))//" and "//TRIM(NUMBER_TO_VSTRING(nfid,"*",ERR,ERROR))&
             & //" fields are not in the same region"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
          ENDIF       
       ENDDO 
     
       !checking whether the field's handle 
       DO nfid=1, FIELDS%NUMBER_OF_FIELDS
          IF(.NOT.ASSOCIATED(FIELDS%FIELDS(nfid)%PTR)) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(nfid,"*",ERR,ERROR))//" field handle in fields list is invalid"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
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
    ENDIF !ASSOCIATED(FIELDS%REGION
   
    !release the pointers if they are associated
    !Allocatable is the refer to dynamic size of array
    IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET)) THEN
       IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)) THEN
          DO nele=1, LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS
             DO ncomp=1, LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%NUMBER_OF_COMPONENTS
                NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS(ncomp)%PTR)
             ENDDO   
             DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET(nele)%COMPONENTS)
          ENDDO
          DEALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET%ELEMENTAL_INFO_SET)
       ENDIF
    ELSE
       ALLOCATE(LOCAL_PROCESS_ELEMENTAL_INFO_SET,STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("could not allocate nodal information set for IO",ERR,ERROR,*999)             
    END IF
   
    !set to number of nodes to zero
    LOCAL_PROCESS_ELEMENTAL_INFO_SET%NUMBER_OF_ELEMENTS=0 
    NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    
    !associated nodel info set with the list of fields
    IF(ASSOCIATED(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)) THEN
       NULLIFY(LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS)
    ENDIF
    LOCAL_PROCESS_ELEMENTAL_INFO_SET%FIELDS=>FIELDS
 
    CALL EXITS("ELEMENTAL_INFO_SET_INITIALISE")
    RETURN
999 CALL ERRORS("ELEMENTAL_INFO_SET_INITIALISE",ERR,ERROR)
    CALL EXITS("ELEMENTAL_INFO_SET_INITIALISE")
    RETURN 1  
  END SUBROUTINE ELEMENTAL_INFO_SET_INITIALISE

  !
  !================================================================================================================================
  !

  !>Export nodal information into multiple files \see{FIELD_IO::EXNODE_INTO_MULTIPLE_FILES}.                 
  SUBROUTINE EXNODE_INTO_MULTIPLE_FILES(FIELDS, FILE_NAME, ERR,ERROR,*)
    !Argument variables       
    TYPE(FIELDS_TYPE), POINTER :: FIELDS !<the field object
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<file name
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    !TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN) :: computational_node_numbers   !<total process number      
    TYPE(FIELD_IO_NODAL_INFO_SET), ALLOCATABLE :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process

    CALL ENTERS("EXNODE_INTO_MULTIPLE_FILES", ERR,ERROR,*999)    

    !Get the number of computational nodes
    computational_node_numbers=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !Get my computational node number
    my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999

    IF(STRCMP(METHOD, "FORTRAIN") THEN
       ALLOCATE(LOCAL_LOCAL_PROCESS_NODAL_INFO_SET),STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodal buffer in IO writing",ERR,ERROR,*999)

       CALL NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FIELDS, ERR,ERROR,*999) 
       CALL NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*999)
       CALL NODAL_INFO_SET_SORT(LOCAL_PROCESS_ELEMENTAL_INFO_SET, my_computational_node_number, ERR,ERROR,*999)    
       CALL EXPORT_NODES_INTO_LOCAL_FILE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, FILE_NAME, my_computational_node_number, computational_node_numbers, &
            &ERR, ERROR, *999)
       CALL NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_ELEMENTAL_INFO_SET, ERR,ERROR,*999)
    ELSE IF(STRCMP(METHOD, "MPIIO") THEN
       CALL FLAG_ERROR("MPIIO is not implemented yet",ERR,ERROR,*999)
    ENDIF         
 
    CALL EXITS("EXNODE_INTO_MULTIPLE_FILES")
    RETURN
999 CALL ERRORS("EXNODE_INTO_MULTIPLE_FILES",ERR,ERROR)
    CALL EXITS("EXNODE_INTO_MULTIPLE_FILES")
    RETURN 1  
  END SUBROUTINE EXNODE_INTO_MULTIPLE_FILES

  !
  !================================================================================================================================
  !  

  !>Initialize nodal information set 
  SUBROUTINE NODAL_INFO_SET_INITIALISE(LOCAL_PROCESS_NODAL_INFO_SET, FIELDS, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), POINTER :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    TYPE(FIELDS_TYPE), POINTER ::FIELDS !<the field object
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: num_field, nn, ncomp  !temporary variable
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("NODAL_INFO_SET_INITIALISE",ERR,ERROR,*999)    

    !validate the input data
    !checking whether the list of fields in the same region
    IF(ASSOCIATED(FIELDS%REGION)) THEN
       DO num_field =2, FIELDS%NUMBER_OF_FIELDS
          IF(FIELDS%FIELDS(num_field-1)%PTR%REGION%USER_NUMBER/=FIELDS%FIELDS(num_field)%PTR%REGION%USER_NUMBER) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(num_field-1,"*",ERR,ERROR))//" and "//TRIM(NUMBER_TO_VSTRING(num_field,"*",ERR,ERROR))&
             & //" fields are not in the same region"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
          ENDIF       
       ENDDO 
       
       !checking whether the field's handle 
       DO num_field=1, FIELDS%NUMBER_OF_FIELDS
          IF(.NOT.ASSOCIATED(FIELDS%FIELDS(num_field)%PTR)) THEN
             LOCAL_ERROR ="No. "//TRIM(NUMBER_TO_VSTRING(num_field,"*",ERR,ERROR))//" field handle in fields list is invalid"
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)         
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
    ENDIF !ASSOCIATED(FIELDS%REGION
   
    !release the pointers if they are associated
    !Allocatable is the refer to dynamic size of array
    IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET)) THEN
       IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
          DO nn=1, LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
             DO ncomp=1, LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS
               NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(ncomp)%PTR)
             ENDDO   
             DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
          ENDDO
          DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)
       ENDIF
    ELSE
       ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET,STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("could not allocate nodal information set for IO",ERR,ERROR,*999)             
    END IF
   
    !set to number of nodes to zero
    LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=0 
    NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    
    !associated nodel info set with the list of fields
    IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)) THEN
      NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)
    ENDIF
    LOCAL_PROCESS_NODAL_INFO_SET%FIELDS=>FIELDS
 
    CALL EXITS("NODAL_INFO_SET_INITIALISE")
    RETURN
999 CALL ERRORS("NODAL_INFO_SET_INITIALISE",ERR,ERROR)
    CALL EXITS("NODAL_INFO_SET_INITIALISE")
    RETURN 1  
  END SUBROUTINE NODAL_INFO_SET_INITIALISE

  !
  !================================================================================================================================
  !  

  !>Collect nodal information from each MPI process      
  SUBROUTINE NODAL_INFO_SET_ATTACH_LOCAL_PROCESS(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), POINTER:: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT):: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx, var_idx, component_idx, np, nn!temporary variable
    LOGICAL :: SWITCH
    INTEGER(INTG), ALLOCATABLE:: NEW_LIST(:)
    TYPE(FIELD_TYPE), POINTER :: FIELD
    !TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !loca domain
    TYPE(DOMAIN_NODES_TYPE), POINTER:: DOMAIN_NODES !nodes in local domain
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER:: tmp_comp(:) !component of field variable
    TYPE(FIELD_VARIABLE_TYPE), POINTER:: FIELD_VARIABLE !field variable

    CALL ENTERS("NODAL_INFO_SET_CREATE_FINISH",ERR,ERROR,*999)    
   
    !initialize the pointer
    NULLIFY(NEW_LIST)
   
    !attache local process to local nodal information set. In current opencmiss system,
    !each local process owns it local nodal information, so all we need to do is to fill the nodal 
    !information set with nodal information of local process 
    IF((LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0).AND.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS) &
    & .AND.(.NOT.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET))) THEN  
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
                            ALLOCATE(NEW_LIST(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                            !add one more node
                            NEW_LIST(1:LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES)=&
                            &LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1:LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES)              
                            NEW_LIST(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1)=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
                            !release the old memory space
                            DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
                            LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>NEW_LIST
                            NULLIFY(NEW_LIST)
                            LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES+1
                         ELSE
                            ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary buffer in IO",ERR,ERROR,*999)
                            LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(1)=DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER
                            LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=1
                         ENDIF   !LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES/=0)                 
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
          NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)                                  
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
                         tmp_comp=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS
                         NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
                      
                         ALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(&
                          &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS+1),STAT=ERR)
                         IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component buffer in IO",ERR,ERROR,*999)
                
                         LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(1:LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)&
                         &%NUMBER_OF_COMPONENTS)=tmp_comp(:)

                         LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)&
                         &%NUMBER_OF_COMPONENTS+1)%PTR=>FIELD_VARIABLE%COMPONENTS(component_idx)                                       
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
       NULLIFY(tmp_comp)
     
       !LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER=>LIST_OF_GLOBAL_NUMBER
       !NULLIFY(LIST_OF_GLOBAL_NUMBER)      
       !release the temporary meomery           
       IF(ASSOCIATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    ELSE
     CALL FLAG_ERROR("nodal information set is not initialized properly, and call start method first",ERR,ERROR,*999)
    ENDIF !(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0).AND.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)
 
    CALL EXITS("NODAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN
999 CALL ERRORS("NODAL_INFO_SET_ATTACH_LOCAL_PROCESS",ERR,ERROR)
    CALL EXITS("NODAL_INFO_SET_ATTACH_LOCAL_PROCESS")
    RETURN 1  
  END SUBROUTINE NODAL_INFO_SET_ATTACH_LOCAL_PROCESS      

  !
  !================================================================================================================================
  !  

  !>Sort nodal information according to the type of field variable component           
  SUBROUTINE NODAL_INFO_SET_SORT(LOCAL_PROCESS_NODAL_INFO_SET, my_computational_node_number, ERR,ERROR,*)      
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), POINTER :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: domain_idx,domain_no, MY_DOMAIN_INDEX !temporary variable
    INTEGER(INTG) :: global_number1, local_number1, global_number2, local_number2
    INTEGER(INTG) :: component_idx, tmp1, nn1, nn2 ! nn, tmp2!temporary variable
    INTEGER(INTG), ALLOCATABLE:: array1(:), array2(:)  
    LOGICAL :: SWITCH
    !INTEGER(INTG), ALLOCATABLE:: LIST_OF_GLOBAL_NUMBER(:)
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER :: tmp_components(:)        
    !TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: tmp_ptr !temporary variable component    
    !TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !field variable
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES1, DOMAIN_NODES2! domain nodes

    !from now on, global numbering are used
    CALL ENTERS("NODAL_INFO_SET_SORT",ERR,ERROR,*999)    
  
    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("list of global numbering in the input data is invalid",ERR,ERROR,*999)     
    ENDIF
    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
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
             tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn2)%COMPONENTS=>&
             &LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%COMPONENTS
             LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn1+1)%COMPONENTS=>tmp_components
                          
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

    CALL EXITS("NODAL_INFO_SET_SORT")
    RETURN
999 CALL ERRORS("NODAL_INFO_SET_SORT",ERR,ERROR)
    CALL EXITS("NODAL_INFO_SET_SORT")
    RETURN 1  
  END SUBROUTINE NODAL_INFO_SET_SORT

  !
  !================================================================================================================================
  !  

  !>Finalize nodal information set    
  SUBROUTINE NODAL_INFO_SET_FINALIZE(LOCAL_PROCESS_NODAL_INFO_SET, ERR,ERROR,*)
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), POINTER :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn, nm  !temporary variable

    CALL ENTERS("NODAL_INFO_SET_FINALIZE",ERR,ERROR,*999)    
   
    IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)) THEN
       NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS)
    ENDIF
    DO nn=1,LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
       IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)) THEN
          DO nm=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS
             IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(nm)%PTR)) THEN
                NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS(nm1)%PTR)
             ENDIF
          ENDDO
          LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%NUMBER_OF_COMPONENTS=0
          LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER=.FALSE.
          IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)) THEN
             DEALLOCATE(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS)
          ENDIF   
       ENDIF
    ENDDO
   
    LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES=0
    IF(ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       NULLIFY(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)
    ENDIF
            
    CALL EXITS("NODAL_INFO_SET_FINALIZE")
    RETURN
999 CALL ERRORS("NODAL_INFO_SET_FINALIZE",ERR,ERROR)
    CALL EXITS("NODAL_INFO_SET_FINALIZE")
    RETURN 1  
  END SUBROUTINE NODAL_INFO_SET_FINALIZE

  !
  !================================================================================================================================
  !  

  !>Get the derivative information               
  FUNCTION LABEL_DERIVATIVE_INFO_GET(GROUP_DERIVATIVES, NUMBER_DERIVATIVES, LABEL_TYPE, ERR, ERROR)
    !Argument variables   
    INTEGER(INTG), INTENT(IN) :: NUMBER_DERIVATIVES
    INTEGER(INTG), INTENT(IN) :: GROUP_DERIVATIVES(NUMBER_DERIVATIVES)
    INTEGER(INTG), INTENT(IN) :: LABEL_TYPE !<identitor for information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dev_idx
    TYPE(VARYING_STRING) ::LABEL_DERIVATIVE_INFO_GET
   
    CALL ENTERS("LABEL_DERIVATIVE_INFO_GET",ERR,ERROR,*999)    
   
    IF(NUMBER_DERIVATIVES==0) THEN
       CALL FLAG_ERROR("number of derivatives in the input data is zero",ERR,ERROR,*999)           
    ENDIF
    IF(LABEL_TYPE/=DERIVATIVE_LABEL) THEN
       CALL FLAG_ERROR("label type in the input data is not derivative label",ERR,ERROR,*999)           
    ENDIF

    IF((NUMBER_DERIVATIVES==1).AND.GROUP_DERIVATIVES(1)==NO_PART_DERIV) THEN
       LABEL_DERIVATIVE_INFO_GET=""
    ELSE
       LABEL_DERIVATIVE_INFO_GET="("
       DO dev_idx=1,NUMBER_DERIVATIVES
          SELECT CASE(GROUP_DERIVATIVES(dev_idx))    
            CASE(NO_PART_DERIV)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET
            CASE(PART_DERIV_S1)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", du/ds1"      
            CASE(PART_DERIV_S1_S1)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds1ds1"       
            CASE(PART_DERIV_S2)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", du/ds2"
            CASE(PART_DERIV_S2_S2)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds2ds2"
            CASE(PART_DERIV_S1_S2)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", du/ds3"
            CASE(PART_DERIV_S3)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds3ds3"
            CASE(PART_DERIV_S3_S3)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds3ds3"
            CASE(PART_DERIV_S1_S3)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds1ds3"
            CASE(PART_DERIV_S2_S3)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds2ds3"
            CASE(PART_DERIV_S1_S2_S3)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^3u/ds1ds2ds3"
            CASE(PART_DERIV_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", du/ds4"
            CASE(PART_DERIV_S4_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds4ds4"
            CASE(PART_DERIV_S1_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds1ds4"
            CASE(PART_DERIV_S2_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds2ds4"
            CASE(PART_DERIV_S3_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^2u/ds3ds4"
            CASE(PART_DERIV_S1_S2_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^3u/ds1ds2ds4"
            CASE(PART_DERIV_S1_S3_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^3u/ds1ds3ds4"
            CASE(PART_DERIV_S2_S3_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^3u/ds2ds3ds4"
            CASE(PART_DERIV_S1_S2_S3_S4)
              LABEL_DERIVATIVE_INFO_GET=LABEL_DERIVATIVE_INFO_GET//", d^4u/ds1ds2ds3ds4"
            CASE DEFAULT
              LABEL_DERIVATIVE_INFO_GET="unknown field variable type, add more details later, #Components="!&
              !&//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
         END SELECT
       ENDDO! dev_idx
    ENDIF !NUMBER_DERIVATIVES==1.AND.GROUP_DERIVATIVES(1)==NO_PART_DERIV
             
    CALL EXITS("LABEL_DERIVATIVE_INFO_GET")
    RETURN
999 CALL ERRORS("LABEL_DERIVATIVE_INFO_GET",ERR,ERROR)
    CALL EXITS("LABEL_DERIVATIVE_INFO_GET")
  END FUNCTION LABEL_DERIVATIVE_INFO_GET    

  !
  !================================================================================================================================
  !  

  !>Get the field information 
  FUNCTION LABEL_FIELD_INFO_GET(COMPONENT, LABEL_TYPE, ERR, ERROR)
    !Argument variables   
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: COMPONENT
    INTEGER(INTG), INTENT(IN) :: LABEL_TYPE !<identitor for information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    TYPE(VARYING_STRING) :: LABEL_FIELD_INFO_GET
   
    CALL ENTERS("LABEL_FIELD_INFO_GET",ERR,ERROR,*999)    
   
    IF(.NOT.ASSOCIATED(COMPONENT)) THEN
       CALL FLAG_ERROR("component pointer in the input data is invalid",ERR,ERROR,*999)           
    ENDIF
   
    IF(LABEL_TYPE/=FIELD_LABEL.AND.LABEL_TYPE/=VARIABLE_LABEL.AND.LABEL_TYPE/=COMPONENT_LABEL) THEN
       CALL FLAG_ERROR("indicator of type in input data is invalid",ERR,ERROR,*999)           
    ENDIF
   
    FIELD=>COMPONENT%FIELD
    VARIABLE=>COMPONENT%FIELD_VARIABLE
    
    SELECT CASE(FIELD%TYPE)
      CASE(FIELD_GEOMETRIC_TYPE) !FIELD_GEOMETRIC_TYPE
        IF(LABEL_TYPE==FIELD_LABEL) THEN
           LABEL_FIELD_INFO_GET="field geometric type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
              !coordinate system              
              SELECT CASE (COMPONENT%REGION%COORDINATE_SYSTEM%TYPE)
                  CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
                   IF(LABEL_TYPE==VARIABLE_LABEL) THEN                
                       LABEL_FIELD_INFO_GET="coordinates,  coordinate, rectangular cartesian"
                   ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                      IF(COMPONENT%COMPONENT_NUMBER==1) THEN
                         LABEL_FIELD_INFO_GET="x"
                      ELSE IF(COMPONENT%COMPONENT_NUMBER==2) THEN
                         LABEL_FIELD_INFO_GET="y"
                      ELSE IF(COMPONENT%COMPONENT_NUMBER==3) THEN
                         LABEL_FIELD_INFO_GET="z"
                  ENDIF                             
                   ENDIF                        
                  !CASE(COORDINATE_CYCLINDRICAL_POLAR_TYPE)
                  !CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
                  !CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
                  !CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
                  CASE DEFAULT
                   IF(LABEL_TYPE==VARIABLE_LABEL) THEN                
                       LABEL_FIELD_INFO_GET="unknown" !coordinates, coordinate, rectangular cartesian,
                   ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                      LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
                   ENDIF                        
              END SELECT
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="Normal_derivative,  field,  normal derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="first_time_derivative,  field,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="second_time_derivative,  field,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown_geometry,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_LABEL   
      CASE(FIELD_FIBRE_TYPE)          
        IF(LABEL_TYPE==FIELD_LABEL) THEN
           LABEL_FIELD_INFO_GET="field fibres type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="fiber,  standand variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="norm_der_fiber,  normal derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="first_time_fiber,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="second_time_fiber,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown_fiber,  unknown field variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_LABEL           
      CASE(FIELD_GENERAL_TYPE)   
        IF(LABEL_TYPE==FIELD_LABEL) THEN
           LABEL_FIELD_INFO_GET="field general type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="general_variabe,  field,  string"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="norm_dev_variable,  field,  string"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="first_time_variable,  field,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="second_time_variable,  field,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown_general,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_LABEL           
      CASE(FIELD_MATERIAL_TYPE)      
        IF(LABEL_TYPE==FIELD_LABEL) THEN
           LABEL_FIELD_INFO_GET="field material type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="material,  field,  standand variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="normal_material,  field,  normal derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="fist_time_material,  field,  firt time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="second_time_material,  field,  second time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown material,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_LABEL           
      CASE DEFAULT
        IF(LABEL_TYPE==FIELD_LABEL) THEN
           LABEL_FIELD_INFO_GET="unknown field type"
        ELSE   
          SELECT CASE(VARIABLE%VARIABLE_TYPE)
            CASE(FIELD_STANDARD_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown,  field,  unknown standand variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_NORMAL_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown,  field,  unknown normal derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV1_VARIABLE_TYPE)
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown,  field,  unknown firt time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE(FIELD_TIME_DERIV2_VARIABLE_TYPE)      
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown, field,  unknown second time derivative of variable"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
            CASE DEFAULT
             IF(LABEL_TYPE==VARIABLE_LABEL) THEN
                LABEL_FIELD_INFO_GET="unknown,  field,  unknown field variable type"
             ELSE IF (LABEL_TYPE==COMPONENT_LABEL) THEN
                LABEL_FIELD_INFO_GET=TRIM(NUMBER_TO_VSTRING(COMPONENT%COMPONENT_NUMBER,"*",ERR,ERROR))
             ENDIF      
          END SELECT!CASE(VARIABLE%VARIABLE_TYPE)
        ENDIF !TYPE==FIELD_LABEL           
    END SELECT
             
    CALL EXITS("LABEL_FIELD_INFO_GET")
    RETURN
999 CALL ERRORS("LABEL_FIELD_INFO_GET",ERR,ERROR)
    CALL EXITS("LABEL_FIELD_INFO_GET")
  END FUNCTION LABEL_FIELD_INFO_GET  

  !
  !================================================================================================================================
  !  

  !>Write the header of a group elements   
  SUBROUTINE EXPORT_NODAL_GROUP_HEADER(LOCAL_PROCESS_NODAL_INFO_SET, LOCAL_NODAL_NUMBER, MAX_NUM_OF_NODAL_DERIVATIVES, &
  &my_computational_node_number, FILE_ID, ERR,ERROR, *)
    !Argument variables  
    TYPE(FIELD_IO_NODAL_INFO_SET), POINTER :: LOCAL_PROCESS_NODAL_INFO_SET  !<LOCAL_PROCESS_NODAL_INFO_SET           
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODAL_NUMBER !<LOCAL_NUMBER IN THE NODAL IO LIST
    INTEGER(INTG), INTENT(INOUT) :: MAX_NUM_OF_NODAL_DERIVATIVES !<MAX_NUM_OF_NODAL_DERIVATIVES
    INTEGER(INTG), INTENT(IN) :: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN) :: FILE_ID !< FILE ID    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUM_OF_FIELDS, NUM_OF_VARIABLES, NUM_OF_NODAL_DEV
    INTEGER(INTG) :: domain_idx, domain_no, MY_DOMAIN_INDEX, local_number, global_number
    INTEGER(INTG), POINTER :: GROUP_FIELDS(:), GROUP_VARIABLES(:), GROUP_DERIVATIVES(:)
    INTEGER(INTG) :: field_idx, comp_idx, comp_idx1, value_idx, var_idx, global_var_idx !dev_idx,   
    TYPE(FIELD_TYPE), POINTER :: field_ptr
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable_ptr   
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER :: tmp_components(:)        
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
    TYPE(VARYING_STRING) :: LINE, LABEL
   
    CALL ENTERS("EXPORT_NODAL_GROUP_HEADER",ERR,ERROR,*999)    
    
    !colllect nodal header information for IO first
    
    !!get the number of this computational node from mpi pool
    !my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)     
    !IF(ERR/=0) GOTO 999        
    
    !attach the temporary pointer
    tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%COMPONENTS    
    
    !collect maximum number of nodal derivatives, number of fields and variables 
    NUM_OF_FIELDS=0
    NUM_OF_VARIABLES=0
    MAX_NUM_OF_NODAL_DERIVATIVES=0
    global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(LOCAL_NODAL_NUMBER)
    NULLIFY(field_ptr)     
    NULLIFY(variable_ptr)
    DO comp_idx=1,LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(LOCAL_NODAL_NUMBER)%NUMBER_OF_COMPONENTS    
       !calculate the number of fields
       IF (.NOT.ASSOCIATED(field_ptr, target=tmp_components(comp_idx)%PTR%FIELD)) THEN
             NUM_OF_FIELDS=NUM_OF_FIELDS+1
           field_ptr=>tmp_components(comp_idx)%PTR%FIELD
       ENDIF         
       
       !calculate the number of variables
       IF (.NOT.ASSOCIATED(variable_ptr, target=tmp_components(comp_idx)%PTR%FIELD_VARIABLE)) THEN
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          variable_ptr=>tmp_components(comp_idx)%PTR%FIELD_VARIABLE
       ENDIF      
   
      !finding the local numbering through the global to local mapping   
       DOMAIN_MAPPING_NODES=>tmp_components(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
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
       IF((.NOT.ASSOCIATED(field_ptr,TARGET=tmp_components(comp_idx)%PTR%FIELD)).AND. &
       &(.NOT.ASSOCIATED(variable_ptr,TARGET=tmp_components(comp_idx)%PTR%FIELD_VARIABLE))) THEN !different field and variables            
          !add one new variable
          NUM_OF_FIELDS=NUM_OF_FIELDS+1
          GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1 
          !add one new component
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
          field_ptr=>tmp_components(comp_idx)%PTR%FIELD
          variable_ptr=>tmp_components(comp_idx)%PTR%FIELD_VARIABLE   
       ELSE IF (ASSOCIATED(field_ptr,TARGET=tmp_components(comp_idx)%PTR%FIELD).AND.&
          &(.NOT.ASSOCIATED(variable_ptr,TARGET=tmp_components(comp_idx)%PTR%FIELD_VARIABLE))) THEN !the same field and  different variables
          !add one new variable
          GROUP_FIELDS(NUM_OF_FIELDS)=GROUP_FIELDS(NUM_OF_FIELDS)+1
          !add one new component
          NUM_OF_VARIABLES=NUM_OF_VARIABLES+1
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1
          variable_ptr=>tmp_components(comp_idx)%PTR%FIELD_VARIABLE   
       ELSE  !different components of the same variable
          !add one new component
          GROUP_VARIABLES(NUM_OF_VARIABLES)=GROUP_VARIABLES(NUM_OF_VARIABLES)+1                               
       ENDIF !field_ptr/=tmp_components%COMPONENTS(comp_idx)%PTR%FIELD     
    ENDDO  !comp_idx        
    
    !write out the nodal header
    var_idx=1
    comp_idx=1
    field_idx=1
    value_idx=1
    comp_idx1=1
    global_var_idx=0
    LINE=" "//"#Fields="//TRIM(NUMBER_TO_VSTRING(SUM(GROUP_FIELDS(1:NUM_OF_FIELDS)),"*",ERR,ERROR))
    CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    DO field_idx=1, NUM_OF_FIELDS              
       !write out the field information
       !LABEL=LABEL_FIELD_INFO_GET(tmp_components(comp_idx1)%PTR, FIELD_LABEL,ERR,ERROR)
       !IF(ERR/=0) THEN
       !   CALL FLAG_ERROR("can not get field label",ERR,ERROR,*999)     
       !   GOTO 999               
       !ENDIF        
       !LINE=TRIM(NUMBER_TO_VSTRING(field_idx,"*",ERR,ERROR))//") "//TRIM(LABEL)&
       !&//" , #variables="//TRIM(NUMBER_TO_VSTRING(GROUP_FIELDS(field_idx),"*",ERR,ERROR))       
       !CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
       
       DO var_idx=1, GROUP_FIELDS(field_idx)
          global_var_idx=global_var_idx+1
          !write out the field information
          LABEL="  "//TRIM(NUMBER_TO_VSTRING(global_var_idx,"*",ERR,ERROR))//") "&
          &//LABEL_FIELD_INFO_GET(tmp_components(comp_idx1)%PTR, VARIABLE_LABEL,ERR,ERROR)
          IF(ERR/=0) THEN
             CALL FLAG_ERROR("can not get variable label",ERR,ERROR,*999)     
             GOTO 999               
          ENDIF        
          LINE=TRIM(LABEL)//", #Components="//TRIM(NUMBER_TO_VSTRING(GROUP_VARIABLES(global_var_idx),"*",ERR,ERROR))        
          CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)              
          
          DO comp_idx=1, GROUP_VARIABLES(global_var_idx)
             !write out the component information
             LABEL="   "//LABEL_FIELD_INFO_GET(tmp_components(comp_idx1)%PTR, COMPONENT_LABEL,ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get component label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF        
             LINE=TRIM(LABEL)//"."                          
             
            !finding the local numbering through the global to local mapping   
             DOMAIN_MAPPING_NODES=>tmp_components(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
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
             LABEL=LABEL_DERIVATIVE_INFO_GET(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV), NUM_OF_NODAL_DEV, DERIVATIVE_LABEL,ERR,ERROR)
             IF(ERR/=0) THEN
                CALL FLAG_ERROR("can not get derivative label",ERR,ERROR,*999)     
                GOTO 999               
             ENDIF
             !assemble the header        
             LINE=LINE//"  Value index= "//TRIM(NUMBER_TO_VSTRING(value_idx,"*",ERR,ERROR))&
             &//", #Derivatives= "//TRIM(NUMBER_TO_VSTRING(NUM_OF_NODAL_DEV-1,"*",ERR,ERROR))//TRIM(LABEL)
             !write out the header             
             CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)
             !increase the component index              
             comp_idx1=comp_idx1+1
             !increase the value index
             value_idx=value_idx+NUM_OF_NODAL_DEV
          ENDDO !comp_idx
       ENDDO !var_idx
    ENDDO !field_idx

    !release temporary memory
    IF(ASSOCIATED(GROUP_FIELDS)) DEALLOCATE(GROUP_FIELDS)
    IF(ASSOCIATED(GROUP_VARIABLES)) DEALLOCATE(GROUP_VARIABLES)
    IF(ASSOCIATED(GROUP_DERIVATIVES)) DEALLOCATE(GROUP_DERIVATIVES)    
        
    CALL EXITS("EXPORT_NODAL_GROUP_HEADER")
    RETURN
999 CALL ERRORS("EXPORT_NODAL_GROUP_HEADER",ERR,ERROR)
    CALL EXITS("EXPORT_NODAL_GROUP_HEADER")
    RETURN 1       
  END SUBROUTINE EXPORT_NODAL_GROUP_HEADER  
 
  !
  !================================================================================================================================
  !  

  !>Write all the nodal information from LOCAL_PROCESS_NODAL_INFO_SET to exnode files       
  SUBROUTINE EXPORT_NODES_INTO_LOCAL_FILE(LOCAL_PROCESS_NODAL_INFO_SET, NAME, my_computational_node_number, computational_node_numbers,ERR, ERROR, *)
    !the reason that my_computational_node_number is used in the argument is for future extension
    !Argument variables   
    TYPE(FIELD_IO_NODAL_INFO_SET), POINTER :: LOCAL_PROCESS_NODAL_INFO_SET !<nodal information in this process
    TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
    INTEGER(INTG), INTENT(IN):: my_computational_node_number !<local process number    
    INTEGER(INTG), INTENT(IN):: computational_node_numbers !<total process number       
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LINE 
    TYPE(VARYING_STRING) :: FILE_NAME !the prefix name of file.    
    INTEGER(INTG) :: FILE_ID, domain_idx, domain_no, local_number, global_number, MY_DOMAIN_INDEX
    INTEGER(INTG), ALLOCATABLE :: GROUP_DERIVATIVES(:)    
    INTEGER(INTG) :: nn, comp_idx,  dev_idx, NUM_OF_NODAL_DEV, MAX_NUM_OF_NODAL_DERIVATIVES
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING_NODES !The domain mapping to calculate nodal mappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES ! domain nodes
    TYPE(FIELD_VARIABLE_COMPONENT_PTR_TYPE), POINTER :: tmp_components(:)            
    REAL(DP), ALLOCATABLE :: NODAL_BUFFER(:), GEOMETRIC_PARAMETERS(:)

    CALL ENTERS("EXPORT_NODES_INTO_LOCAL_FILE",ERR,ERROR,*999)    
    
    !get my own computianal node number--be careful the rank of process in the MPI pool 
    !is not necessarily equal to numbering of computional node, so use method COMPUTATIONAL_NODE_NUMBER_GET
    !will be a secured way to get the number         
    FILE_NAME=NAME//".part"//TRIM(NUMBER_TO_VSTRING(my_computational_node_number,"*",ERR,ERROR))//".exnode"        
    FILE_ID=1029
    MAX_NUM_OF_NODAL_DERIVATIVES=0
    
    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET)) THEN
       CALL FLAG_ERROR("the nodal information set in input is invalid",ERR,ERROR,*999)            
    ENDIF

    IF(.NOT.ASSOCIATED(LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER)) THEN
       CALL FLAG_ERROR("the nodal information set is not associated with any numbering list",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES==0) THEN 
       CALL FLAG_ERROR("the nodal information set does not contain any nodes",ERR,ERROR,*999)            
    ENDIF

    IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(1)%SAME_HEADER==.TRUE.) THEN
       CALL FLAG_ERROR("the first header flag of nodal information set should be false",ERR,ERROR,*999)            
    ENDIF
    
    NULLIFY(NODAL_BUFFER)    
    NULLIFY(GROUP_DERIVATIVES)
    
    !open a file 
    CALL FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, ERR,ERROR,*999)
    
    !write out the group name    
    LINE=FILEDS_GROUP_INFO_GET(LOCAL_PROCESS_NODAL_INFO_SET%FIELDS, ERR,ERROR)    
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get group namein IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF               
    CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       
    !write out the number of files    
    LINE=MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)      
    IF(ERR/=0) THEN
       CALL FLAG_ERROR("can not get multiple file information in IO",ERR,ERROR,*999)     
       GOTO 999               
    ENDIF             
    !CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)       

    DO nn=1, LOCAL_PROCESS_NODAL_INFO_SET%NUMBER_OF_NODES
      
       tmp_components=>LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%COMPONENTS         
       global_number=LOCAL_PROCESS_NODAL_INFO_SET%LIST_OF_GLOBAL_NUMBER(nn)
       
       !check whether need to write out the nodal information header  
       IF(LOCAL_PROCESS_NODAL_INFO_SET%NODAL_INFO_SET(nn)%SAME_HEADER==.FALSE.) THEN
        !write out the nodal header
           
          CALL EXPORT_NODAL_GROUP_HEADER(LOCAL_PROCESS_NODAL_INFO_SET, nn, MAX_NUM_OF_NODAL_DERIVATIVES, &
          &my_computational_node_number, FILE_ID, ERR,ERROR,*999) 
           !value_idx=value_idx-1 !the len of NODAL_BUFFER                      
          !checking: whether need to allocate temporary memory for Io writing
          IF(ASSOCIATED(NODAL_BUFFER)) THEN
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
          DOMAIN_MAPPING_NODES=>tmp_components(comp_idx)%PTR%DOMAIN%MAPPINGS%NODES       
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
             LINE="Node:     "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NODES%NODES(local_number)%USER_NUMBER,"*",ERR,ERROR))
             CALL FORTRAN_FILE_WRITE_STRING(FILE_ID, LINE, LEN_TRIM(LINE), ERR,ERROR,*999)                      
          ENDIF
          
          !get the nodal partial derivatives
          NUM_OF_NODAL_DEV=DOMAIN_NODES%NODES(local_number)%NUMBER_OF_DERIVATIVES
          GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV)=DOMAIN_NODES%NODES(local_number)%PARTIAL_DERIVATIVE_INDEX(:)               
          !sort  the partial derivatives
          CALL LIST_SORT(GROUP_DERIVATIVES(1:NUM_OF_NODAL_DEV),ERR,ERROR,*999)
          DO dev_idx=1, NUM_OF_NODAL_DEV
             NULLIFY(GEOMETRIC_PARAMETERS)
             IF(comp_idx==1) THEN
                CALL FIELD_PARAMETER_SET_GET(tmp_components(comp_idx)%PTR%FIELD,&
                                          &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
             ELSE IF (ASSOCIATED(tmp_components(comp_idx)%PTR%DOMAIN, TARGET=tmp_components(comp_idx-1)%PTR%DOMAIN)) THEN
                CALL FIELD_PARAMETER_SET_GET(tmp_components(comp_idx)%PTR%FIELD,&
                                          &FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
             ENDIF                                       
           NODAL_BUFFER(dev_idx)=GEOMETRIC_PARAMETERS(tmp_components(comp_idx)%PTR%&
            &PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(GROUP_DERIVATIVES(dev_idx),local_number,0))
               CALL FORTRAN_FILE_WRITE_DP(FILE_ID, NODAL_BUFFER, NUM_OF_NODAL_DEV, ERR,ERROR,*999)                
          ENDDO
       ENDDO !comp_idx                   
       
    ENDDO!nn    
    
    !close a file    
    CALL FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*999)

    !release the temporary memory
    IF(ASSOCIATED(NODAL_BUFFER)) THEN
       DEALLOCATE(NODAL_BUFFER)
    ENDIF
    IF(ASSOCIATED(GROUP_DERIVATIVES)) THEN
       DEALLOCATE(GROUP_DERIVATIVES)
    ENDIF
    
    CALL EXITS("EXPORT_NODES_INTO_LOCAL_FILE")
    RETURN
999 CALL ERRORS("EXPORT_NODES_INTO_LOCAL_FILE",ERR,ERROR)
    CALL EXITS("EXPORT_NODES_INTO_LOCAL_FILE")
    RETURN 1  
  END SUBROUTINE EXPORT_NODES_INTO_LOCAL_FILE
    
  !
  !================================================================================================================================
  !  

  !>Get the region label  
  FUNCTION FILEDS_GROUP_INFO_GET(FIELDS, ERR, ERROR)
    !Argument variables   
    TYPE(FIELDS_TYPE), INTENT(IN):: FIELDS !<FILEDS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: FILEDS_GROUP_INFO_GET !<On exit group name
   
    CALL ENTERS("FILEDS_GROUP_INFO_GET",ERR,ERROR,*999)    
   
    IF(ASSOCIATED(FIELDS%REGION)) THEN
       FILEDS_GROUP_INFO_GET="Group name: "//TRIM(FIELDS%REGION%LABEL)
    ELSE
       CALL FLAG_ERROR("fields are not associated with a region, so can not find a group name",ERR,ERROR,*999)     
    ENDIF                
     
    CALL EXITS("FILEDS_GROUP_INFO_GET")
    RETURN
999 CALL ERRORS("FILEDS_GROUP_INFO_GET",ERR,ERROR)
    CALL EXITS("FILEDS_GROUP_INFO_GET")
  END FUNCTION FILEDS_GROUP_INFO_GET

  !
  !================================================================================================================================
  !  

  !>Get the number of files      
  FUNCTION MULTI_FILES_INFO_GET(computational_node_numbers, ERR, ERROR)
    !Argument variables   
    INTEGER(INTG), INTENT(IN) :: computational_node_numbers  
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   !Temporary variables
    TYPE(VARYING_STRING) :: MULTI_FILES_INFO_GET !<ON exit multi fiels information
   
   CALL ENTERS("MULTI_FILES_INFO_GET",ERR,ERROR,*999)    
   
    MULTI_FILES_INFO_GET="files: "//TRIM(NUMBER_TO_VSTRING(computational_node_numbers,"*",ERR,ERROR)); 
     
    CALL EXITS("MULTI_FILES_INFO_GET")
    RETURN
999 CALL ERRORS("MULTI_FILES_INFO_GET",ERR,ERROR)
    CALL EXITS("MULTI_FILES_INFO_GET")
  END FUNCTION MULTI_FILES_INFO_GET  

  !
  !================================================================================================================================
  !  

  !>Write a string     
  SUBROUTINE FORTRAN_FILE_WRITE_STRING(FILE_ID, STRING_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(IN) :: STRING_DATA !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    IF(LEN_OF_DATA==0) THEN
       CALL FLAG_ERROR("leng of string is zero",ERR,ERROR,*999) 
    ENDIF
    CALL ENTERS("FORTRAN_FILE_WRITE_STRING",ERR,ERROR,*999)    
   
    WRITE(FILE_ID, "(A)"), CHAR(STRING_DATA) 
    CALL EXITS("FORTRAN_FILE_WRITE_STRING")
    RETURN
999 CALL ERRORS("FORTRAN_FILE_WRITE_STRING",ERR,ERROR)
    CALL EXITS("FORTRAN_FILE_WRITE_STRING")
    RETURN 1
  END SUBROUTINE FORTRAN_FILE_WRITE_STRING

  !
  !================================================================================================================================
  !  

  !>Write a real data      
  SUBROUTINE FORTRAN_FILE_WRITE_DP(FILE_ID, REAL_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    REAL(DP), INTENT(IN) :: REAL_DATA(:) !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: DP_FMT !<the name of file.    
    
    CALL ENTERS("FORTRAN_FILE_WRITE_DP",ERR,ERROR,*999)    
   
    DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"ES)"
    WRITE(FILE_ID, CHAR(DP_FMT)), REAL_DATA(1:LEN_OF_DATA)
    
    CALL EXITS("FORTRAN_FILE_WRITE_DP")
    RETURN
999 CALL ERRORS("FORTRAN_FILE_WRITE_DP",ERR,ERROR)
    CALL EXITS("FORTRAN_FILE_WRITE_DP")
    RETURN 1
  END SUBROUTINE FORTRAN_FILE_WRITE_DP

  !
  !================================================================================================================================
  !  

  !>Write a integer data      
  SUBROUTINE FORTRAN_FILE_WRITE_INTG(FILE_ID, INTG_DATA, LEN_OF_DATA, ERR,ERROR,*)
  
    !Argument variables   
    INTEGER(INTG), INTENT(IN) :: INTG_DATA(:) !<the name of file.
    INTEGER(INTG), INTENT(IN) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(IN) :: LEN_OF_DATA !<length of string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: DP_FMT !<the name of file. 
       
   CALL ENTERS("FORTRAN_FILE_WRITE_INTG",ERR,ERROR,*999)    
   
   DP_FMT="("//TRIM(NUMBER_TO_VSTRING(LEN_OF_DATA,"*",ERR,ERROR))//"I)"
   WRITE(FILE_ID, CHAR(DP_FMT)), INTG_DATA(1:LEN_OF_DATA)
    
    CALL EXITS("FORTRAN_FILE_WRITE_INTG")
    RETURN
999 CALL ERRORS("FORTRAN_FILE_WRITE_INTG",ERR,ERROR)
    CALL EXITS("FORTRAN_FILE_WRITE_INTG")
    RETURN 1
  END SUBROUTINE FORTRAN_FILE_WRITE_INTG
 
  !
  !================================================================================================================================
  !  

  !>Open a file using Fortran 
  SUBROUTINE FORTRAN_FILE_OPEN(FILE_ID, FILE_NAME, ERR,ERROR,*)
  
    !Argument variables   
    TYPE(VARYING_STRING), INTENT(INOUT) :: FILE_NAME !<the name of file.
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
        
   CALL ENTERS("FORTRAN_FILE_OPEN",ERR,ERROR,*999)       

   OPEN(UNIT=FILE_ID, FILE=CHAR(FILE_NAME), STATUS='REPLACE', FORM="FORMATTED", ERR=999)   
    
    CALL EXITS("FORTRAN_FILE_OPEN")
    RETURN
999 CALL ERRORS("FORTRAN_FILE_OPEN",ERR,ERROR)
    CALL EXITS("FORTRAN_FILE_OPEN")
    RETURN 1
  END SUBROUTINE FORTRAN_FILE_OPEN
  !
  !================================================================================================================================
  !  

  !>Close a file using Fortran 
  SUBROUTINE FORTRAN_FILE_CLOSE(FILE_ID, ERR,ERROR,*)
  
    !Argument variables   
    INTEGER(INTG), INTENT(INOUT) :: FILE_ID !<file ID
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
   CALL ENTERS("FORTRAN_FILE_CLOSE",ERR,ERROR,*999)    
   
   CLOSE(UNIT=FILE_ID, ERR=999)   
    
    CALL EXITS("FORTRAN_FILE_CLOSE")
    RETURN
999 CALL ERRORS("FORTRAN_FILE_CLOSE",ERR,ERROR)
    CALL EXITS("FORTRAN_FILE_CLOSE")
    RETURN 1
  END SUBROUTINE FORTRAN_FILE_CLOSE
 
END MODULE FIELD_IO_ROUTINES

   
