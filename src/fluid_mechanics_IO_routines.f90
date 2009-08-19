MODULE FLUID_MECHANICS_IO_ROUTINES

   USE BASE_ROUTINES
!   USE LISTS
!   USE BASIS_ROUTINES
!   USE MESH_ROUTINES
!   USE NODE_ROUTINES
!   USE COMP_ENVIRONMENT
!   USE COORDINATE_ROUTINES
!   USE ISO_VARYING_STRING
!   USE REGION_ROUTINES
!   USE MACHINE_CONSTANTS
!   USE KINDS
   USE FIELD_ROUTINES
!   USE ISO_VARYING_STRING
!   USE ISO_C_BINDING
!   USE STRINGS
   USE TYPES
!   USE CONSTANTS
!   USE MPI
!   USE CMISS_MPI
!   USE INPUT_OUTPUT

  IMPLICIT NONE

  !1=M, 2=V, 3=P !


  INTEGER, DIMENSION(:), ALLOCATABLE:: NodesPerElement
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::ElementNodesScales
  INTEGER, DIMENSION(:,:), ALLOCATABLE::ElementNodes
  INTEGER:: NumberOfFields
  INTEGER:: NumberOfDimensions
  INTEGER:: ValueIndex
  INTEGER:: NumberOfVariableComponents
  INTEGER:: NumberOfMeshComponents
  INTEGER:: NumberOfMaterialComponents
  INTEGER:: NumberOfNodesDefined
  INTEGER:: NumberOfFieldComponent(3)
  INTEGER:: NumberOfElements
  INTEGER:: GlobalElementNumber(10)
  INTEGER:: MaxNodesPerElement
  INTEGER:: MaxNodesPerMeshComponent

  INTEGER:: ELEMENT_NUMBER
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: XI_COORDINATES,COORDINATES
   DOUBLE PRECISION:: test

  INTEGER:: lagrange_simplex

  INTEGER, DIMENSION(:), ALLOCATABLE:: NodesPerMeshComponent

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeXValue
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeYValue
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeZValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeUValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeVValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeWValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodePValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeMUValue
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeRHOValue  

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeUValue_analytic
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeVValue_analytic 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeWValue_analytic 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodePValue_analytic 

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeUValue_error
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeVValue_error
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeWValue_error
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodePValue_error

  INTEGER, DIMENSION(:), ALLOCATABLE::SimplexOutputHelp


  TYPE(FIELD_TYPE), POINTER :: FIELD
  TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
  TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT

  DOUBLE PRECISION:: ScaleFactorsPerElementNodes(10,10)

  DOUBLE PRECISION:: MaxNodeUValue_error, MaxNodeVValue_error, MaxNodeWValue_error, MaxNodePValue_error

  INTEGER:: TRI_BASIS, TET_BASIS, QUAD_BASIS, HEX_BASIS
  CHARACTER*2 NMs(99),KNOT


  !1=M, 2=V, 3=P !
        TYPE ARRAY_MESH
          INTEGER Lx                                    ! External Structure: Number of Tessellation Coefficients
          INTEGER Lt                                    ! External Structure: Number of Tessellation elements
          INTEGER ID                                    ! Internal Structure: Basis number to call letter map
          INTEGER Lb
          DOUBLE PRECISION, pointer :: X(:,:)           ! External Structure: Tessellation coefficients
          INTEGER, pointer :: T(:,:)                    ! External Structure: Tessellation map
          INTEGER, pointer :: NBT(:,:)                  ! Internal Structure: Neighboring element map
          INTEGER, pointer :: B(:,:)
        END TYPE ARRAY_MESH

        TYPE ARRAY_INT
          DOUBLE PRECISION, pointer :: Y(:,:)           ! Internal Structure: Basis @ given volume quadrature points
          DOUBLE PRECISION, pointer :: Y_f(:,:,:)       ! Internal Structure: Basis @ given facet quadrature points
          DOUBLE PRECISION, pointer :: dY(:,:,:), TdY(:,:,:), dY_f(:,:,:,:)
        END TYPE ARRAY_INT

        TYPE ARRAY_BASE
          INTEGER n                                     ! External Structure: Function Space dimension on T_e
          INTEGER nl                                    ! Internal Structure: Breakdown for tensor input
          INTEGER DISCONT                               ! External Structure: Continuity constraint
          INTEGER DM                                    ! External Structure: Field dimension
          TYPE(ARRAY_INT) I                             ! External Structure: Field integrations
          DOUBLE PRECISION, pointer :: Q(:,:)           ! External Structure: Basis coefficient tensor
          DOUBLE PRECISION, pointer :: P(:,:)           ! External Structure: Basis power tensor
          DOUBLE PRECISION, pointer :: XI(:,:)          ! External Structure: Basis xi-point coordinates
          INTEGER, pointer :: B_ID(:,:)                 ! External Structure: Basis ordering tensor
          CHARACTER*2 CL                                ! External Structure: Basis call letter
        END TYPE ARRAY_BASE

        TYPE ARRAY_PROBLEM_BASE
          TYPE(ARRAY_BASE), pointer :: B(:)             ! Internal Structure: Basis
          INTEGER n_pts                                 ! External Structure: number of volume quadrature points
          INTEGER n_pts_f                               ! External Structure: number of facet quadrature points
          INTEGER n_ptsl                                ! Internal Structure: volume quadrature for tensor input
          INTEGER n_pts_fl                              ! Internal Structure: facet quadrature for tensor input
          INTEGER HEXA                                  ! External Structure: tensor input parameter
          INTEGER FACES                                 ! External Structure: number of master element facets
          INTEGER FNODES                                ! External Structure: number of nodes per facet (spatial map)
          INTEGER n_B                                   ! External Structure: number of basis
          INTEGER DM                                    ! External Structure: spatial dimension of master element
          INTEGER SPL, VEL, PRS
          INTEGER TRI_BASIS, QUAD_BASIS, TET_BASIS, HEX_BASIS
          DOUBLE PRECISION VL                           ! External Structure: volume of master element
          DOUBLE PRECISION VL_f(6)                      ! External Structure: facet volume of master element
          DOUBLE PRECISION nrm(3,6)                     ! Internal Structure: facet normals of master element
          DOUBLE PRECISION, pointer :: gpt(:,:)         ! External Structure: volume quadrature points
          DOUBLE PRECISION, pointer :: gw(:)            ! External Structure: volume quadrature weights
          DOUBLE PRECISION, pointer :: gpt_f(:,:,:)     ! External Structure: facet quadrature points
          DOUBLE PRECISION, pointer :: gw_f(:)          ! External Structure: facet quadrature weights
        END TYPE ARRAY_PROBLEM_BASE

       TYPE EXPORT_CONTAINER
               DOUBLE PRECISION, POINTER::N(:,:)
               INTEGER, POINTER::M(:,:),V(:,:),P(:,:)
               INTEGER:: D,F,ID_M,ID_V,ID_P,IT_M,IT_V,IT_P,IT_T
               INTEGER:: E_M,E_P,E_V,E_T,EN_M,EN_P,EN_V,EN_T,N_M,N_P,N_V,N_T
       END TYPE EXPORT_CONTAINER


       TYPE DARCY_PARAMETERS
               INTEGER:: TESTCASE
               INTEGER:: BC_NUMBER_OF_WALL_NODES, NUMBER_OF_BCS
               DOUBLE PRECISION:: PERM, VIS, PERM_OVER_VIS
               DOUBLE PRECISION:: X1, X2, Y1, Y2, Z1, Z2
               DOUBLE PRECISION:: LENGTH, GEOM_TOL
               DOUBLE PRECISION:: max_node_spacing
               LOGICAL :: STAB, DEBUG, ANALYTIC
       END TYPE DARCY_PARAMETERS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  TYPE (ARRAY_PROBLEM_BASE) BASE_INFO
  TYPE (ARRAY_MESH) MESH_INFO(3)
  TYPE (DARCY_PARAMETERS) DARCY
  INTEGER FLD, DIMEN, OPENCMISS_INTERPOLATION(3),a,b
  INTEGER NumberOfNodesPerElement(3), ArrayOfNodesDefined(3), NumberOfElementsDefined(3), TotalNumberOfNodes
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::OPENCMISS_NODE_COORD
  INTEGER, DIMENSION(:,:), ALLOCATABLE::OPENCMISS_ELEM_M,OPENCMISS_ELEM_V,OPENCMISS_ELEM_P
  CHARACTER*60 IN_CHAR
  CHARACTER*90 NIMZ
!   CHARACTER*30 NAMz
  CHARACTER*90 NAMz
  INTEGER:: ALLOC_ERROR


  CONTAINS

! --------------------------------------------------------------------------------------------------------------------
! 
! chrm, 19.08.2009
! 
  SUBROUTINE READ_DARCY_PARAMS

  IMPLICIT NONE

  CHARACTER*90 DARCY_PARAM_FILE

  DARCY_PARAM_FILE='./input/Darcy_parameters.inp'

  WRITE(*,*)'Reading Darcy parameters.'

  OPEN(UNIT=37,FILE=DARCY_PARAM_FILE,STATUS='old',action='read') ! Read base file for initial parameters


      DO WHILE (0 < 1)
        READ(37,*,END=50) IN_CHAR
        IF (INDEX(IN_CHAR,'TESTCASE:') == 1)      READ(37,*) DARCY%TESTCASE

        IF (INDEX(IN_CHAR,'STAB:') == 1)          READ(37,*) DARCY%STAB
        IF (INDEX(IN_CHAR,'ANALYTIC:') == 1)      READ(37,*) DARCY%ANALYTIC
        IF (INDEX(IN_CHAR,'DEBUG:') == 1)         READ(37,*) DARCY%DEBUG

        IF (INDEX(IN_CHAR,'LENGTH:') == 1)        READ(37,*) DARCY%LENGTH

        IF (INDEX(IN_CHAR,'GEOM_TOL:') == 1)      READ(37,*) DARCY%GEOM_TOL
        IF (INDEX(IN_CHAR,'X1:') == 1)            READ(37,*) DARCY%X1
        IF (INDEX(IN_CHAR,'X2:') == 1)            READ(37,*) DARCY%X2
        IF (INDEX(IN_CHAR,'Y1:') == 1)            READ(37,*) DARCY%Y1
        IF (INDEX(IN_CHAR,'Y2:') == 1)            READ(37,*) DARCY%Y2
        IF (INDEX(IN_CHAR,'Z1:') == 1)            READ(37,*) DARCY%Z1
        IF (INDEX(IN_CHAR,'Z2:') == 1)            READ(37,*) DARCY%Z2

        IF (INDEX(IN_CHAR,'PERM:') == 1)          READ(37,*) DARCY%PERM
        IF (INDEX(IN_CHAR,'VIS:') == 1)           READ(37,*) DARCY%VIS

        IF (INDEX(IN_CHAR,'BC_NUMBER_OF_WALL_NODES:') == 1) READ(37,*) DARCY%BC_NUMBER_OF_WALL_NODES
        IF (INDEX(IN_CHAR,'NUMBER_OF_BCS:') == 1) READ(37,*) DARCY%NUMBER_OF_BCS
      END DO

50 CLOSE(37)

  IF( DARCY%DEBUG ) THEN
    write(*,*)'Read Darcy parameters from the following file = ',DARCY_PARAM_FILE
    write(*,*)'Press ENTER to continue.'
    read(*,*)
  END IF

  IF( ABS(DARCY%VIS) > 1.0E-14 ) THEN
    DARCY%PERM_OVER_VIS = DARCY%PERM / DARCY%VIS
  ELSE
    WRITE(*,*)'Darcy_parameters: VIS cannot be machine zero.'
    STOP
  END IF

  DARCY%max_node_spacing = 2.0

  END SUBROUTINE READ_DARCY_PARAMS
! 
! --------------------------------------------------------------------------------------------------------------------


  !HERE THE TYPES DEFINED ABOVE ARE FILLED WITH THE DATA PROVIDED  
  SUBROUTINE FLUID_MECHANICS_IO_WRITE_CMGUI(REGION, NAME, ERR, ERROR,*)

  INTEGER:: I,J,K,L,M,N

  INTEGER :: ERR !<The error code
  TYPE(VARYING_STRING):: ERROR !<The error string

   TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for
   TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.

   IF (ALLOCATED(NodesPerElement)) DEALLOCATE(NodesPerElement)
   IF (ALLOCATED(NodesPerMeshComponent)) DEALLOCATE(NodesPerMeshComponent)
   IF (ALLOCATED(XI_COORDINATES)) DEALLOCATE(XI_COORDINATES)
   IF (ALLOCATED(COORDINATES)) DEALLOCATE(COORDINATES)
   IF (ALLOCATED(NodeXValue)) DEALLOCATE(NodeXValue)
   IF (ALLOCATED(NodeYValue)) DEALLOCATE(NodeYValue)
   IF (ALLOCATED(NodeZValue)) DEALLOCATE(NodeZValue)
   IF (ALLOCATED(NodeUValue)) DEALLOCATE(NodeUValue)
   IF (ALLOCATED(NodeVValue)) DEALLOCATE(NodeVValue)
   IF (ALLOCATED(NodeWValue)) DEALLOCATE(NodeWValue)
   IF (ALLOCATED(NodePValue)) DEALLOCATE(NodePValue)
   IF (ALLOCATED(NodeMUValue)) DEALLOCATE(NodeMUValue)
   IF (ALLOCATED(NodeRHOValue)) DEALLOCATE(NodeRHOValue)
   IF (ALLOCATED(ElementNodesScales)) DEALLOCATE(ElementNodesScales)
   IF (ALLOCATED(ElementNodes)) DEALLOCATE(ElementNodes)

       KNOT = '0'
       NMs(1) = '1'
       NMs(2) = '2'
       NMs(3) = '3'
       NMs(4) = '4'
       NMs(5) = '5'
       NMs(6) = '6'
       NMs(7) = '7'
       NMs(8) = '8'
       NMs(9) = '9'
       K = 9
       DO I = 1,9
          K = K + 1
          NMs(K) = TRIM(NMs(I))//TRIM(KNOT)
          DO J = 1,9
             K = K + 1
             NMs(K) = TRIM(NMs(I))//TRIM(NMs(J))
          END DO
       END DO




   NumberOfFields=REGION%fields%number_of_fields
   NumberOfDimensions=REGION%coordinate_system%number_of_dimensions

  NumberOfVariableComponents=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%&
  &variables(1)%number_of_components

   NumberOfMaterialComponents=REGION%equations_sets%equations_sets(1)%ptr%materials%materials_field%&
   &variables(1)%number_of_components

   NumberOfElements=REGION%meshes%meshes(1)%ptr%number_of_elements
   NumberOfMeshComponents=REGION%meshes%meshes(1)%ptr%number_of_components
   ALLOCATE(NodesPerElement(NumberOfMeshComponents))

   ALLOCATE(NodesPerMeshComponent(NumberOfMeshComponents))
   MaxNodesPerElement=0
   DO I=1,NumberOfMeshComponents
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1)&
      &%ptr%topology%elements%elements(1)%basis%number_of_element_parameters
      NodesPerMeshComponent(I)=REGION%meshes%meshes(1)%ptr%topology(I)%ptr%nodes%number_of_nodes
   END DO


   MaxNodesPerElement=NodesPerElement(1)



   MaxNodesPerMeshComponent=NodesPerMeshComponent(1)

   ALLOCATE(XI_COORDINATES(NumberOfDimensions))
   ALLOCATE(COORDINATES(NumberOfDimensions))

   ALLOCATE(NodeXValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeYValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeZValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeUValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeVValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeWValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodePValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeMUValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeRHOValue(NodesPerMeshComponent(1)))
   ALLOCATE(ElementNodesScales(NumberOfElements,NodesPerElement(1)))
   ALLOCATE(ElementNodes(NumberOfElements,NodesPerElement(1)))

   ALLOCATE(NodeUValue_analytic(NodesPerMeshComponent(1)))
   ALLOCATE(NodeVValue_analytic(NodesPerMeshComponent(1)))
   ALLOCATE(NodeWValue_analytic(NodesPerMeshComponent(1)))
   ALLOCATE(NodePValue_analytic(NodesPerMeshComponent(1)))

   ALLOCATE(NodeUValue_error(NodesPerMeshComponent(1)))
   ALLOCATE(NodeVValue_error(NodesPerMeshComponent(1)))
   ALLOCATE(NodeWValue_error(NodesPerMeshComponent(1)))
   ALLOCATE(NodePValue_error(NodesPerMeshComponent(1)))

! THIS NEEDS TO BE ADJUSTED NOW!!!!!!



    CALL ENTERS("CMGUI OUTPUT",ERR,ERROR,*999)



   FIELD=>REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field



   NULLIFY(INTERPOLATION_PARAMETERS)
   NULLIFY(INTERPOLATED_POINT)
   CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,FIELD_U_VARIABLE_TYPE,INTERPOLATION_PARAMETERS&
   &,ERR,ERROR,*999)



   CALL FIELD_INTERPOLATED_POINT_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)



  DO I=1,NumberOfElements
   DO J=1,NodesPerElement(1)
 
      ELEMENT_NUMBER=I
      XI_COORDINATES(1)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%&
      &geometric_interp_parameters%bases(1)%ptr%node_position_index(J,1)-1.0)/(REGION%equations_sets%&
      &equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1)&
      &%ptr%number_of_nodes_xi(1)-1.0)
      XI_COORDINATES(2)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%&
      &geometric_interp_parameters%bases(1)%ptr%node_position_index(J,2)-1.0)/(REGION%equations_sets%&
      &equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1)&
      &%ptr%number_of_nodes_xi(2)-1.0)
      IF(NumberOfDimensions==3)THEN
      XI_COORDINATES(3)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%&
      &geometric_interp_parameters%bases(1)%ptr%node_position_index(J,3)-1.0)/(REGION%equations_sets%&
      &equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1)&
      &%ptr%number_of_nodes_xi(3)-1.0)
      END IF



      !K is global node number
      K=REGION%meshes%meshes(1)%ptr%topology(1)%ptr%elements%elements(I)%global_element_nodes(J)



      IF(NumberOfDimensions==3)THEN
      COORDINATES=(/1,1,1/)
      ELSE IF(NumberOfDimensions==2)THEN
      COORDINATES=(/1,1/)
      END IF



      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,&
      &INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
 



      CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,INTERPOLATED_POINT,ERR,ERROR,*999)

      NodeXValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)




      NodeYValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))

      IF(NumberOfDimensions==3)THEN
      NodeZValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
      END IF



      NodeUValue(K)=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)
      NodeVValue(K)=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))
      IF(NumberOfDimensions==3)THEN
      NodeWValue(K)=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))
      END IF

! ! !       NodeUValue(K)=INTERPOLATED_POINT%VALUES(1,1)
! ! !       NodeVValue(K)=INTERPOLATED_POINT%VALUES(2,1)
! ! !       NodeWValue(K)=INTERPOLATED_POINT%VALUES(3,1)

      IF(NumberOfDimensions==3)THEN
      NodePValue(K)=INTERPOLATED_POINT%VALUES(4,1)
      ELSE IF(NumberOfDimensions==2)THEN
      NodePValue(K)=INTERPOLATED_POINT%VALUES(3,1)
      END IF



   END DO 
  END DO




   NodeMUValue=REGION%equations_sets%equations_sets(1)%ptr%materials%materials_field%variables(1)%&
   &parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1)

   NodeRHOValue=REGION%equations_sets%equations_sets(1)%ptr%materials%materials_field%variables(1)%&
   &parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(2)




   IF( NumberOfDimensions==3 )THEN
     !For 3D, the following call works ...
     lagrange_simplex=REGION%equations_sets%equations_sets(1)%ptr%equations%&
     &interpolation%geometric_interp_parameters%bases(1)%ptr%type
   ELSE
     ! ... but the above call does not work for 2D.
     !Thus, for 2D, we hard-wire it to 'quad':
     lagrange_simplex=1
   END IF


                   

  NumberOfFieldComponent(1)=NumberOfDimensions
  NumberOfFieldComponent(2)=NumberOfVariableComponents
  NumberOfFieldComponent(3)=NumberOfMaterialComponents




    DO I=1,NumberOfElements
    DO J=1,NodesPerElement(1)
      ElementNodes(I,J)=REGION%meshes%meshes(1)%ptr%topology(1)%&
      &ptr%elements%elements(I)%global_element_nodes(J)
      ElementNodesScales(I,J)=1.0000000000000000E+00
    END DO
    END DO




     CALL WRITE_NODES_CMGUI(NAME)

     CALL WRITE_ELEMENTS_CMGUI(NAME)



    CALL EXITS("CMGUI OUTPUT")
    RETURN
999 CALL ERRORS("CMGUI OUTPUT",ERR,ERROR)


  END SUBROUTINE FLUID_MECHANICS_IO_WRITE_CMGUI


! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

  

  SUBROUTINE WRITE_NODES_CMGUI(NAME)

  IMPLICIT NONE

  TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
  TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.


  INTEGER:: I,J,K,L,M,N



  DOUBLE PRECISION:: COORD_X, COORD_Y, COORD_Z, ARG_X, ARG_Y, ARG_Z
  DOUBLE PRECISION:: FACT

  IF( DARCY%TESTCASE == 1 ) THEN
    FACT = DARCY%PERM_OVER_VIS
  ELSE
    FACT = 2.0_DP * PI * DARCY%PERM_OVER_VIS / DARCY%LENGTH
!     FACT = 1.0_DP
  END IF


  IF( DARCY%ANALYTIC ) THEN
    IF( NumberOfDimensions==2 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        COORD_X = NodeXValue(I)
        COORD_Y = NodeYValue(I)
        ARG_X = 2.0_DP * PI * COORD_X / DARCY%LENGTH
        ARG_Y = 2.0_DP * PI * COORD_Y / DARCY%LENGTH

        IF( DARCY%TESTCASE == 1 ) THEN
          NodeUValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X + 2.0_DP * COORD_Y )
          NodeVValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X - 2.0_DP * COORD_Y )
          NodePValue_analytic(I) = COORD_X * COORD_X + 2.0_DP * COORD_X * COORD_Y - COORD_Y * COORD_Y
        ELSE
          NodeUValue_analytic(I) = - FACT * COS( ARG_X ) * SIN( ARG_Y ) 
          NodeVValue_analytic(I) = - FACT * SIN( ARG_X ) * COS( ARG_Y ) 
          NodePValue_analytic(I) =          SIN( ARG_X ) * SIN( ARG_Y )
        END IF

        NodeUValue_error(I) = NodeUValue(I) - NodeUValue_analytic(I)
        NodeVValue_error(I) = NodeVValue(I) - NodeVValue_analytic(I)
        NodePValue_error(I) = NodePValue(I) - NodePValue_analytic(I)
      END DO
    ELSE IF( NumberOfDimensions==3 ) THEN
      DO I = 1,NodesPerMeshComponent(1)

        COORD_X = NodeXValue(I)
        COORD_Y = NodeYValue(I)
        COORD_Z = NodeZValue(I)

        ARG_X = 2.0_DP * PI * COORD_X / DARCY%LENGTH
        ARG_Y = 2.0_DP * PI * COORD_Y / DARCY%LENGTH
        ARG_Z = 2.0_DP * PI * COORD_Z / DARCY%LENGTH

        IF( DARCY%TESTCASE == 1 ) THEN
          NodeUValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X + 2.0_DP * COORD_Y + COORD_Z )
          NodeVValue_analytic(I) = - FACT * ( 2.0_DP * COORD_X - 2.0_DP * COORD_Y + COORD_Z )
          NodeWValue_analytic(I) = - FACT * ( 3.0_DP + COORD_X + COORD_Y )
          NodePValue_analytic(I) = COORD_X * COORD_X + 2.0_DP * COORD_X * COORD_Y - COORD_Y * COORD_Y + &
            & 3.0_DP * COORD_Z + COORD_Z * COORD_X  + COORD_Z * COORD_Y 
        ELSE
          NodeUValue_analytic(I) = - FACT * COS( ARG_X ) * SIN( ARG_Y )  * SIN( ARG_Z ) 
          NodeVValue_analytic(I) = - FACT * SIN( ARG_X ) * COS( ARG_Y )  * SIN( ARG_Z )  
          NodeWValue_analytic(I) = - FACT * SIN( ARG_X ) * SIN( ARG_Y )  * COS( ARG_Z )  
          NodePValue_analytic(I) =          SIN( ARG_X ) * SIN( ARG_Y )  * SIN( ARG_Z )  
        END IF

        NodeUValue_error(I) = NodeUValue(I) - NodeUValue_analytic(I)
        NodeVValue_error(I) = NodeVValue(I) - NodeVValue_analytic(I)
        NodeWValue_error(I) = NodeWValue(I) - NodeWValue_analytic(I)
        NodePValue_error(I) = NodePValue(I) - NodePValue_analytic(I)

      END DO

    END IF

  END IF



       FILENAME="./output/"//NAME//".exnode"

       OPEN(UNIT=14, FILE=CHAR(FILENAME),STATUS='unknown')




! WRITING HEADER INFORMATION

       WRITE(14,*) 'Group name: OpenCMISS'
       IF( DARCY%ANALYTIC ) THEN
         WRITE(14,*) '#Fields=',TRIM(NMs(NumberOfFields + 2))
       ELSE
         WRITE(14,*) '#Fields=',TRIM(NMs(NumberOfFields))
       END IF
       ValueIndex=1
       WRITE(14,*) ' 1) coordinates,  coordinate, rectangular cartesian, #Components=',TRIM(NMs(NumberOfDimensions))
       DO I=1,NumberOfDimensions
         IF(I==1) THEN
           WRITE(14,*) '   x.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
         ELSE IF(I==2) THEN
           WRITE(14,*) '   y.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
         ELSE
           WRITE(14,*) '   z.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
         END IF
         ValueIndex=ValueIndex+1
       END DO
       WRITE(14,*) ' 2) general,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))
       DO I=1,NumberOfVariableComponents
         WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
         ValueIndex=ValueIndex+1
       END DO
       WRITE(14,*) ' 3) material,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfMaterialComponents))
       DO I=1,NumberOfMaterialComponents
         WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
         ValueIndex=ValueIndex+1
       END DO

       IF( DARCY%ANALYTIC ) THEN
         WRITE(14,*) ' 4) exact,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))
         DO I=1,NumberOfVariableComponents
           WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
           ValueIndex=ValueIndex+1
         END DO

         WRITE(14,*) ' 5) error,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))
         DO I=1,NumberOfVariableComponents
           WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
           ValueIndex=ValueIndex+1
         END DO
       END IF

! NOW WRITE NODE INFORMATION
! ! ! 
            DO I = 1,NodesPerMeshComponent(1)
               WRITE(14,*) ' Node: ',I

                  WRITE(14,'("    ", es25.16 )')NodeXValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeYValue(I)
                  IF(NumberOfDimensions==3) THEN
                  WRITE(14,'("    ", es25.16 )')NodeZValue(I)
                  END IF
                  WRITE(14,'("    ", es25.16 )')NodeUValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeVValue(I)
                  IF(NumberOfDimensions==3) THEN
                  WRITE(14,'("    ", es25.16 )')NodeWValue(I)
                  END IF
                  WRITE(14,'("    ", es25.16 )')NodePValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeMUValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeRHOValue(I)

                  IF( DARCY%ANALYTIC ) THEN
                    WRITE(14,'("    ", es25.16 )')NodeUValue_analytic(I)
                    WRITE(14,'("    ", es25.16 )')NodeVValue_analytic(I)
                    IF(NumberOfDimensions==3) THEN
                      WRITE(14,'("    ", es25.16 )')NodeWValue_analytic(I)
                    END IF
                    WRITE(14,'("    ", es25.16 )')NodePValue_analytic(I)

                    WRITE(14,'("    ", es25.16 )')NodeUValue_error(I)
                    WRITE(14,'("    ", es25.16 )')NodeVValue_error(I)
                    IF(NumberOfDimensions==3) THEN
                      WRITE(14,'("    ", es25.16 )')NodeWValue_error(I)
                    END IF
                    WRITE(14,'("    ", es25.16 )')NodePValue_error(I)
                  END IF
            END DO


       WRITE(14,*) ' '
       CLOSE(14)

  WRITE(*,*)'Writing Nodes...'

  ! ----------------------------------------------------------
  ! Write file to monitor convergence of discretization error
  IF( DARCY%ANALYTIC ) THEN
    OPEN(UNIT=23, FILE='./output/conv.node',STATUS='unknown')

    MaxNodeUValue_error = 0.0
    MaxNodeVValue_error = 0.0
    MaxNodeWValue_error = 0.0
    MaxNodePValue_error = 0.0

    IF( NumberOfDimensions==2 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        IF( abs(mod( ((NodeXValue(I)-DARCY%X1) / DARCY%max_node_spacing), 1.0)) < DARCY%GEOM_TOL ) THEN
          IF( abs(mod( ((NodeYValue(I)-DARCY%Y1) / DARCY%max_node_spacing), 1.0)) < DARCY%GEOM_TOL ) THEN

              WRITE(23,'("    ", es25.16 )')NodeXValue(I)
              WRITE(23,'("    ", es25.16 )')NodeYValue(I)

              WRITE(23,'("    ", es25.16 )')NodeUValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodeVValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodePValue_error(I)

              WRITE(23,*) ' '

              IF( abs(NodeUValue_error(I)) > MaxNodeUValue_error ) MaxNodeUValue_error = abs(NodeUValue_error(I))
              IF( abs(NodeVValue_error(I)) > MaxNodeVValue_error ) MaxNodeVValue_error = abs(NodeVValue_error(I))
              IF( abs(NodePValue_error(I)) > MaxNodePValue_error ) MaxNodePValue_error = abs(NodePValue_error(I))

          END IF
        END IF
      END DO
      WRITE(23,'("    MaxNodeUValue_error = ", es25.16 )')MaxNodeUValue_error
      WRITE(23,'("    MaxNodeVValue_error = ", es25.16 )')MaxNodeVValue_error
      WRITE(23,'("    MaxNodePValue_error = ", es25.16 )')MaxNodePValue_error
      WRITE(23,*) ' '
    ELSE IF( NumberOfDimensions==3 ) THEN
      DO I = 1,NodesPerMeshComponent(1)
        IF( abs(mod( ((NodeXValue(I)-DARCY%X1) / DARCY%max_node_spacing), 1.0)) < DARCY%GEOM_TOL ) THEN
          IF( abs(mod( ((NodeYValue(I)-DARCY%Y1) / DARCY%max_node_spacing), 1.0)) < DARCY%GEOM_TOL ) THEN
            IF( abs(mod( ((NodeZValue(I)-DARCY%Z1) / DARCY%max_node_spacing), 1.0)) < DARCY%GEOM_TOL ) THEN

              WRITE(23,'("    ", es25.16 )')NodeXValue(I)
              WRITE(23,'("    ", es25.16 )')NodeYValue(I)
              WRITE(23,'("    ", es25.16 )')NodeZValue(I)

              WRITE(23,'("    ", es25.16 )')NodeUValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodeVValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodeWValue_error(I)
              WRITE(23,'("    ", es25.16 )')NodePValue_error(I)

              WRITE(23,*) ' '

              IF( abs(NodeUValue_error(I)) > MaxNodeUValue_error ) MaxNodeUValue_error = abs(NodeUValue_error(I))
              IF( abs(NodeVValue_error(I)) > MaxNodeVValue_error ) MaxNodeVValue_error = abs(NodeVValue_error(I))
              IF( abs(NodeWValue_error(I)) > MaxNodeWValue_error ) MaxNodeWValue_error = abs(NodeWValue_error(I))
              IF( abs(NodePValue_error(I)) > MaxNodePValue_error ) MaxNodePValue_error = abs(NodePValue_error(I))

            END IF
          END IF
        END IF
      END DO
      WRITE(23,'("    MaxNodeUValue_error = ", es25.16 )')MaxNodeUValue_error
      WRITE(23,'("    MaxNodeVValue_error = ", es25.16 )')MaxNodeVValue_error
      IF( NumberOfDimensions==3 ) THEN
        WRITE(23,'("    MaxNodeWValue_error = ", es25.16 )')MaxNodeWValue_error
      END IF
      WRITE(23,'("    MaxNodePValue_error = ", es25.16 )')MaxNodePValue_error
      WRITE(23,*) ' '
    END IF

    CLOSE(23)
  END IF
  ! ----------------------------------------------------------

  END SUBROUTINE WRITE_NODES_CMGUI


! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

  SUBROUTINE WRITE_ELEMENTS_CMGUI(NAME)

   IMPLICIT NONE


  TYPE(VARYING_STRING), INTENT(IN) :: NAME !<the prefix name of file.
  TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.

  CHARACTER*60 ELEM_TYPE
  INTEGER:: I,J,K,L,M,N

       FILENAME="./output/"//NAME//".exelem"

       OPEN(UNIT=5, FILE=CHAR(FILENAME),STATUS='unknown')

       WRITE(5,*) 'Group name: OpenCMISS'

       IF(lagrange_simplex==2) THEN
         WRITE(5,*) 'Shape.  Dimension=',TRIM(NMs(NumberOfDimensions)),', simplex(2;3)*simplex*simplex'
         IF(MaxNodesPerElement==3) THEN
           WRITE(5,*) '#Scale factor sets= 1'
           WRITE(5,*) ' l.simplex(2)*l.simplex, #Scale factors= ', NodesPerElement(1)
         ELSE IF(MaxNodesPerElement==4) THEN
           WRITE(5,*) '#Scale factor sets= 1'
           WRITE(5,*) ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', NodesPerElement(1)
          ELSE IF (MaxNodesPerElement== 10 ) THEN
           WRITE(5,*) '#Scale factor sets= 1'
           WRITE(5,*) ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', NodesPerElement(1)
         ELSE
           WRITE(5,*) '#Scale factor sets= 0'
         END IF


       ELSE IF (lagrange_simplex==1) THEN
         WRITE(5,*) 'Shape.  Dimension= ',TRIM(NMs(NumberOfDimensions))
         WRITE(5,*) '#Scale factor sets= 1'
         IF(NumberOfDimensions==2) THEN
             IF(MaxNodesPerElement==4) THEN
                   WRITE(5,*) 'l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==9) THEN
                   WRITE(5,*) 'q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==16) THEN
                   WRITE(5,*) 'c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(1)
             END IF
         ELSE
             IF(MaxNodesPerElement==8) THEN
                   WRITE(5,*) 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==27) THEN
                   WRITE(5,*) 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==64) THEN
                   WRITE(5,*) 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(1)
             END IF
         END IF
      END IF

       WRITE(5,*) '#Nodes= ',TRIM(NMs(NodesPerElement(1)))
       WRITE(5,*) '#Fields= ',TRIM(Nms(NumberOfFields))


       DO I=1,NumberOfFields


          IF(I==1)THEN
           WRITE(5,*)' 1) coordinates,  coordinate, rectangular cartesian, #Components= ',TRIM(NMs(NumberOfDimensions))
          ELSE IF(I==2) THEN
           WRITE(5,*)' 2) general,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfVariableComponents))
          ELSE IF(I==3) THEN
           WRITE(5,*)' 3) material,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfMaterialComponents))
          END IF

      DO J=1,NumberOfFieldComponent(I)

        IF(NumberOfDimensions==2) THEN

             IF(I==1)THEN
              IF(J==1) THEN
                IF(MaxNodesPerElement==4)THEN
                    WRITE(5,*)'   x.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9) THEN
                    WRITE(5,*)'   x.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                    WRITE(5,*)'   x.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF 
               ELSE IF(J==2) THEN
                IF(MaxNodesPerElement==4) THEN
                    WRITE(5,*)'   y.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                    WRITE(5,*)'   y.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                    WRITE(5,*)'   y.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
               ELSE IF(J==3) THEN
                IF(MaxNodesPerElement==4) THEN
                    WRITE(5,*)'   z.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                    WRITE(5,*)'   z.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                    WRITE(5,*)'   z.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
              END IF
             ELSE
                IF(MaxNodesPerElement==4) THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
             END IF


          ELSE IF(NumberOfDimensions==3) THEN

             IF(I==1)THEN
              IF(J==1) THEN
                IF(MaxNodesPerElement==8) THEN
                    WRITE(5,*)'   x.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                    WRITE(5,*)'   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                    WRITE(5,*)'   x.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'

                ELSE IF(MaxNodesPerElement==4)  THEN
                    WRITE(5,*)'   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                    WRITE(5,*)'   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                    WRITE(5,*)'   x.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'

                END IF 
               ELSE IF(J==2) THEN
                IF(MaxNodesPerElement==8) THEN
                    WRITE(5,*)'   y.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                    WRITE(5,*)'   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                    WRITE(5,*)'   y.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'

                ELSE IF(MaxNodesPerElement==4)  THEN
                    WRITE(5,*)'   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                    WRITE(5,*)'   y.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                    WRITE(5,*)'   y.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'


                END IF
               ELSE IF(J==3) THEN
                IF(MaxNodesPerElement==8) THEN
                    WRITE(5,*)'   z.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                    WRITE(5,*)'   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                    WRITE(5,*)'   z.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'

                ELSE IF(MaxNodesPerElement==4)  THEN
                    WRITE(5,*)'   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                    WRITE(5,*)'   z.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                    WRITE(5,*)'   z.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'



                END IF
              END IF
             ELSE
                IF(MaxNodesPerElement==8) THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'

                ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'


                END IF
             END IF
          END IF






               WRITE(5,*) '   #Nodes= ',TRIM(NMs(MaxNodesPerElement))



           DO K = 1,MaxNodesPerElement
               WRITE(5,*) '    ',TRIM(NMs(K)),'.  #Values=1'
               WRITE(5,*) '     Value indices:     1'
               WRITE(5,*) '     Scale factor indices:   ',TRIM(NMs(K))
           END DO

          END DO

        END DO

       IF(lagrange_simplex==2) THEN

       ALLOCATE(SimplexOutputHelp(NodesPerElement(1)))

       DO K = 1,NumberOfElements
         IF(NumberOfDimensions==2)THEN
              SimplexOutputHelp=ElementNodes(K,1:NodesPerElement(1))
         ELSE IF(NumberOfDimensions==3) THEN
              SimplexOutputHelp(1)=ElementNodes(K,1)
              SimplexOutputHelp(2)=ElementNodes(K,5)
              SimplexOutputHelp(3)=ElementNodes(K,2)
              SimplexOutputHelp(4)=ElementNodes(K,7)
              SimplexOutputHelp(5)=ElementNodes(K,10)
              SimplexOutputHelp(6)=ElementNodes(K,4)
              SimplexOutputHelp(7)=ElementNodes(K,6)
              SimplexOutputHelp(8)=ElementNodes(K,8)
              SimplexOutputHelp(9)=ElementNodes(K,9)
              SimplexOutputHelp(10)=ElementNodes(K,3)
         END IF
            WRITE(5,*) 'Element:     ', K,' 0  0'
            WRITE(5,*) '   Nodes:'
            WRITE(5,*) '   ', SimplexOutputHelp
            WRITE(5,*) '   Scale factors:'
            WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(1))
        END DO

       ELSE IF (lagrange_simplex==1) THEN

       DO K = 1,NumberOfElements
            WRITE(5,*) 'Element:     ', K,' 0  0'
            WRITE(5,*) '   Nodes:'
            WRITE(5,*) '   ', ElementNodes(K,1:NodesPerElement(1))
            WRITE(5,*) '   Scale factors:'
            WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(1))
       END DO

       END IF

       WRITE(5,*) ' '
       CLOSE(5)



  WRITE(*,*)'Writing Elements...'

  END SUBROUTINE WRITE_ELEMENTS_CMGUI


! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

  SUBROUTINE FLUID_MECHANICS_IO_READ_CMHEART(EXPORT,ERR,ERROR,*)

  IMPLICIT NONE

  TYPE (EXPORT_CONTAINER):: EXPORT  
  TYPE (EXPORT_CONTAINER):: TMP  
  INTEGER :: ERR !<The error code
  TYPE(VARYING_STRING):: ERROR !<The error string
  INTEGER:: test
  INTEGER:: I,J,K,L,M,N

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*)'Run from bin directory if input needed.'
  WRITE(*,*)
  WRITE(*,*)'REMEMBER: M,V,P need to be defined as required by cmHeart!'
  WRITE(*,*)
  WRITE(*,*)'Press ENTER to start'
  READ(*,*)

  OPEN(UNIT=42, FILE='./input/CMHEART.inp',STATUS='old')
!  OPEN(UNIT=42, FILE='./CMHEART.inp',STATUS='old')


  CALL READ_AUX
  CALL READ_NODES
  CALL READ_ELEMENTS


  ALLOCATE(OPENCMISS_ELEM_M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)),STAT=ALLOC_ERROR)
  ALLOCATE(OPENCMISS_ELEM_V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)),STAT=ALLOC_ERROR)
  ALLOCATE(OPENCMISS_ELEM_P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)),STAT=ALLOC_ERROR)


  CALL MAKE_UNIQUE

  CALL ORDER_NUMBERING(OPENCMISS_ELEM_M,MESH_INFO(1)%T,NumberOfElementsDefined(1),NumberOfNodesPerElement(1),1)
  CALL ORDER_NUMBERING(OPENCMISS_ELEM_V,MESH_INFO(2)%T,NumberOfElementsDefined(2),NumberOfNodesPerElement(2),2) 
  CALL ORDER_NUMBERING(OPENCMISS_ELEM_P,MESH_INFO(3)%T,NumberOfElementsDefined(3),NumberOfNodesPerElement(3),3)  
  
 

  ! CALL PRINT_ON_SCREEN  

  WRITE(*,*)'export finished successfully...'
  WRITE(*,*)

  IF(ALLOC_ERROR.NE.0) THEN
    STOP 'Error during allocation'
  END IF



  ALLOCATE(TMP%M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)),STAT=ALLOC_ERROR)
  ALLOCATE(TMP%V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)),STAT=ALLOC_ERROR)
  ALLOCATE(TMP%P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)),STAT=ALLOC_ERROR)
  ALLOCATE(TMP%N(TotalNumberOfNodes,3),STAT=ALLOC_ERROR)


  TMP%M=OPENCMISS_ELEM_M

  TMP%V=OPENCMISS_ELEM_V
  TMP%P=OPENCMISS_ELEM_P

  TMP%N=OPENCMISS_NODE_COORD

  TMP%D=DIMEN
  TMP%F=BASE_INFO%n_B

  TMP%ID_M=1
  TMP%ID_V=2
  TMP%ID_P=3

  TMP%IT_M=OPENCMISS_INTERPOLATION(1)
  TMP%IT_V=OPENCMISS_INTERPOLATION(2)
  TMP%IT_P=OPENCMISS_INTERPOLATION(3)
  
  IF (BASE_INFO%HEXA==1) THEN
    !LAGRANGIAN BASIS
    TMP%IT_T=1
  ELSE 
    ! SIMPLEX BASIS
    TMP%IT_T=2
  END IF

  TMP%E_M=NumberOfElementsDefined(1)
  TMP%E_V=NumberOfElementsDefined(2)
  TMP%E_P=NumberOfElementsDefined(3)
  TMP%E_T=NumberOfElementsDefined(3)

  TMP%EN_M=NumberOfNodesPerElement(1)
  TMP%EN_V=NumberOfNodesPerElement(2)
  TMP%EN_P=NumberOfNodesPerElement(3)
  TMP%EN_T=TMP%EN_M+TMP%EN_V+TMP%EN_P


  TMP%N_M=ArrayOfNodesDefined(1)
  TMP%N_V=ArrayOfNodesDefined(2)
  TMP%N_P=ArrayOfNodesDefined(3)
  TMP%N_T=TMP%N_M+TMP%N_V+TMP%N_P

  EXPORT=TMP

  IF(ALLOC_ERROR.NE.0) THEN
    STOP 'Error during allocation'
  END IF


  END SUBROUTINE FLUID_MECHANICS_IO_READ_CMHEART

! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------


SUBROUTINE READ_AUX
      IMPLICIT NONE

  INTEGER:: I,J,K,L,M,N

  READ(42,*) NIMZ

  NIMZ = TRIM(NIMZ); BASE_INFO%n_B = 0; BASE_INFO%HEXA = 0; BASE_INFO%DM = 3

    OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read') ! Read base file for initial parameters

      DO WHILE (0 < 1)
        READ(1,*,END=50) IN_CHAR
        IF (INDEX(IN_CHAR,'no_fields!') == 1)         READ(1,*) BASE_INFO%n_B
        IF (INDEX(IN_CHAR,'no_gauss!') == 1)          READ(1,*) BASE_INFO%n_pts
        IF (INDEX(IN_CHAR,'volume!') == 1)            READ(1,*) BASE_INFO%VL
        IF (INDEX(IN_CHAR,'no_gauss_f!') == 1)        READ(1,*) BASE_INFO%n_pts_f
        IF (INDEX(IN_CHAR,'no_ele_faces!') == 1)      READ(1,*) BASE_INFO%FACES
        IF (INDEX(IN_CHAR,'no_ele_nodes_f!') == 1)    READ(1,*) BASE_INFO%FNODES
        IF (INDEX(IN_CHAR,'hexa_basis!') == 1)        BASE_INFO%HEXA = 1
        IF (INDEX(IN_CHAR,'domain_dimension!') == 1)  READ(1,*) BASE_INFO%DM
        IF (INDEX(IN_CHAR,'TRI_BASIS!') == 1)  BASE_INFO%TRI_BASIS = 1
        IF (INDEX(IN_CHAR,'TET_BASIS!') == 1)  BASE_INFO%TET_BASIS = 1
        IF (INDEX(IN_CHAR,'QUAD_BASIS!') == 1) BASE_INFO%QUAD_BASIS = 1
        IF (INDEX(IN_CHAR,'HEX_BASIS!') == 1)  BASE_INFO%HEX_BASIS = 1
      END DO

50 CLOSE(1)

    DIMEN=BASE_INFO%DM


! ! ! !         WRITE(*,*)'BASE_INFO%n_B',BASE_INFO%n_B
! ! ! !         WRITE(*,*)'BASE_INFO%n_pts',BASE_INFO%n_pts
! ! ! !         WRITE(*,*)'BASE_INFO%VL',BASE_INFO%VL
! ! ! !         WRITE(*,*)'BASE_INFO%n_pts_f',BASE_INFO%n_pts_f
! ! ! !         WRITE(*,*)'BASE_INFO%FACES',BASE_INFO%FACES
! ! ! !         WRITE(*,*)'BASE_INFO%FNODES',BASE_INFO%FNODES
! ! ! !         WRITE(*,*)'BASE_INFO%HEXA',BASE_INFO%HEXA
! ! ! !         WRITE(*,*)'BASE_INFO%DM',BASE_INFO%DM
! ! ! !         WRITE(*,*)'BASE_INFO%TRI_BASIS',BASE_INFO%TRI_BASIS
! ! ! !         WRITE(*,*)'BASE_INFO%TET_BASIS',BASE_INFO%TET_BASIS
! ! ! !         WRITE(*,*)'BASE_INFO%QUAD_BASIS',BASE_INFO%QUAD_BASIS
! ! ! !         WRITE(*,*)'BASE_INFO%HEX_BASIS',BASE_INFO%HEX_BASIS

        ALLOCATE(BASE_INFO%B(BASE_INFO%n_B),STAT=ALLOC_ERROR)


    OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read') ! Read base file for initial parameters

      DO WHILE (0 < 1)
        READ(1,*,END=52) IN_CHAR
        IF (INDEX(IN_CHAR,'no_basis_M!') == 1) READ(1,*) BASE_INFO%B(1)%n
        IF (INDEX(IN_CHAR,'no_basis_V!') == 1) READ(1,*) BASE_INFO%B(2)%n
        IF (INDEX(IN_CHAR,'no_basis_P!') == 1) READ(1,*) BASE_INFO%B(3)%n

        IF (INDEX(IN_CHAR,'dim_field_M!') == 1) READ(1,*) BASE_INFO%B(1)%DM
        IF (INDEX(IN_CHAR,'dim_field_V!') == 1) READ(1,*) BASE_INFO%B(2)%DM
        IF (INDEX(IN_CHAR,'dim_field_P!') == 1) READ(1,*) BASE_INFO%B(3)%DM
      END DO
52 CLOSE(1)


IF (BASE_INFO%QUAD_BASIS/= 1.AND.BASE_INFO%HEX_BASIS /= 1.AND.BASE_INFO%TRI_BASIS/= 1.AND.BASE_INFO%TET_BASIS /= 1)THEN
  
    WRITE(*,*)'Cubic Hermite not implemented yet'
ELSE

DO I=1,3

        IF (BASE_INFO%TRI_BASIS== 1.OR.BASE_INFO%TET_BASIS == 1)THEN
            IF(DIMEN==2) THEN
              IF(BASE_INFO%B(I)%n==3) THEN
                    NumberOfNodesPerElement(I)=3
                    OPENCMISS_INTERPOLATION(I)=7
              ELSE IF(BASE_INFO%B(I)%n==6) THEN
                    NumberOfNodesPerElement(I)=6
                    OPENCMISS_INTERPOLATION(I)=8
              ELSE IF(BASE_INFO%B(I)%n==10) THEN
                    NumberOfNodesPerElement(I)=10
                    OPENCMISS_INTERPOLATION(I)=9
              ELSE
                  STOP
              END IF
           ELSE IF(DIMEN==3) THEN
              IF(BASE_INFO%B(I)%n==4) THEN
                    NumberOfNodesPerElement(I)=4
                    OPENCMISS_INTERPOLATION(I)=7
              ELSE IF(BASE_INFO%B(I)%n==10) THEN
                    NumberOfNodesPerElement(I)=10
                    OPENCMISS_INTERPOLATION(I)=8
              ELSE IF(BASE_INFO%B(I)%n==20) THEN
                    NumberOfNodesPerElement(I)=20
                    OPENCMISS_INTERPOLATION(I)=9
              ELSE
                  STOP
              END IF
              ELSE 
                STOP
           END IF
        END IF
 

        IF (BASE_INFO%QUAD_BASIS== 1.OR.BASE_INFO%HEX_BASIS == 1)THEN
            IF(BASE_INFO%B(I)%n==2) THEN
            !2D/3D LINEAR LAGRANGE
                  OPENCMISS_INTERPOLATION(I)=1
                  IF(DIMEN==2) THEN
                    NumberOfNodesPerElement(I)=4
                  ELSE IF(DIMEN==3) THEN
                    NumberOfNodesPerElement(I)=8
                  ELSE 
                    STOP
                  END IF



             ELSE IF(BASE_INFO%B(I)%n==3) THEN
            !2D/3D QUADRATIC LAGRANGE
                  OPENCMISS_INTERPOLATION(I)=2
                  IF(DIMEN==2) THEN
                    NumberOfNodesPerElement(I)=9
                  ELSE IF(DIMEN==3) THEN
                    NumberOfNodesPerElement(I)=27
                  ELSE 
                    STOP
                  END IF



            ELSE IF(BASE_INFO%B(I)%n==4) THEN
            !2D/3D CUBIC LAGRANGE
                  OPENCMISS_INTERPOLATION(I)=3
                  IF(DIMEN==2) THEN
                    NumberOfNodesPerElement(I)=16
                  ELSE IF(DIMEN==3) THEN
                    NumberOfNodesPerElement(I)=64
                  ELSE 
                    STOP
                  END IF

            ELSE
                  STOP
            END IF
           END IF

END DO



END IF

  IF(ALLOC_ERROR.NE.0) THEN
    STOP 'Error during allocation'
  END IF

END SUBROUTINE READ_AUX

! ----------------------------------------------------------------------------------

SUBROUTINE ORDER_NUMBERING(NEW,OLD,n,m,I)

      IMPLICIT NONE
  INTEGER:: I,J,K,L,M,N
  INTEGER::NEW(n,m),OLD(n,m)


!  NEW=OLD

      DO J=1,n
      IF (BASE_INFO%QUAD_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==2) THEN
        !2D HEX LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==3) THEN
        !2D HEX QUADR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,5)
        NEW(J,3)=OLD(J,2)
        NEW(J,4)=OLD(J,6)
        NEW(J,5)=OLD(J,7)
        NEW(J,6)=OLD(J,8)
        NEW(J,7)=OLD(J,3)
        NEW(J,8)=OLD(J,9)
        NEW(J,9)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==4) THEN
        !2D HEX CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,5)
        NEW(J,3)=OLD(J,6)
        NEW(J,4)=OLD(J,2)
        NEW(J,5)=OLD(J,7)
        NEW(J,6)=OLD(J,8)
        NEW(J,7)=OLD(J,9)
        NEW(J,8)=OLD(J,10)
        NEW(J,9)=OLD(J,11)
        NEW(J,10)=OLD(J,12)
        NEW(J,11)=OLD(J,13)
        NEW(J,12)=OLD(J,14)
        NEW(J,13)=OLD(J,3)
        NEW(J,14)=OLD(J,15)
        NEW(J,15)=OLD(J,16)
        NEW(J,16)=OLD(J,4)
        ELSE
           STOP
        END IF

      ELSE IF (BASE_INFO%HEX_BASIS == 1) THEN


        IF(BASE_INFO%B(I)%n==2) THEN
        !3D HEX LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,6)
        NEW(J,7)=OLD(J,7)
        NEW(J,8)=OLD(J,8)



        ELSE IF(BASE_INFO%B(I)%n==3) THEN
        !3D HEX QUADR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,9)
        NEW(J,3)=OLD(J,2)
        NEW(J,4)=OLD(J,10)
        NEW(J,5)=OLD(J,11)
        NEW(J,6)=OLD(J,12)
        NEW(J,7)=OLD(J,3)
        NEW(J,8)=OLD(J,13)
        NEW(J,9)=OLD(J,4)
        NEW(J,10)=OLD(J,14)
        NEW(J,11)=OLD(J,15)
        NEW(J,12)=OLD(J,16)
        NEW(J,13)=OLD(J,17)
        NEW(J,14)=OLD(J,18)
        NEW(J,15)=OLD(J,19)
        NEW(J,16)=OLD(J,20)
        NEW(J,17)=OLD(J,21)
        NEW(J,18)=OLD(J,22)
        NEW(J,19)=OLD(J,5)
        NEW(J,20)=OLD(J,23)
        NEW(J,21)=OLD(J,6)
        NEW(J,22)=OLD(J,24)
        NEW(J,23)=OLD(J,25)
        NEW(J,24)=OLD(J,26)
        NEW(J,25)=OLD(J,7)
        NEW(J,26)=OLD(J,27)
        NEW(J,27)=OLD(J,8)


        ELSE IF(BASE_INFO%B(I)%n==4) THEN
    !3D HEX CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,9)
        NEW(J,3)=OLD(J,10)
        NEW(J,4)=OLD(J,2)
        NEW(J,5)=OLD(J,11)
        NEW(J,6)=OLD(J,12)
        NEW(J,7)=OLD(J,13)
        NEW(J,8)=OLD(J,14)
        NEW(J,9)=OLD(J,15)
        NEW(J,10)=OLD(J,16)
        NEW(J,11)=OLD(J,17)
        NEW(J,12)=OLD(J,18)
        NEW(J,13)=OLD(J,3)
        NEW(J,14)=OLD(J,19)
        NEW(J,15)=OLD(J,20)
        NEW(J,16)=OLD(J,4)
        NEW(J,17)=OLD(J,21)
        NEW(J,18)=OLD(J,22)
        NEW(J,19)=OLD(J,23)
        NEW(J,20)=OLD(J,24)
        NEW(J,21)=OLD(J,25)
        NEW(J,22)=OLD(J,26)
        NEW(J,23)=OLD(J,27)
        NEW(J,24)=OLD(J,28)
        NEW(J,25)=OLD(J,29)
        NEW(J,26)=OLD(J,30)
        NEW(J,27)=OLD(J,31)
        NEW(J,28)=OLD(J,32)
        NEW(J,29)=OLD(J,33)
        NEW(J,30)=OLD(J,34)
        NEW(J,31)=OLD(J,35)
        NEW(J,32)=OLD(J,36)
        NEW(J,33)=OLD(J,37)
        NEW(J,34)=OLD(J,38)
        NEW(J,35)=OLD(J,39)
        NEW(J,36)=OLD(J,40)
        NEW(J,37)=OLD(J,41)
        NEW(J,38)=OLD(J,42)
        NEW(J,39)=OLD(J,43)
        NEW(J,40)=OLD(J,44)
        NEW(J,41)=OLD(J,45)
        NEW(J,42)=OLD(J,46)
        NEW(J,43)=OLD(J,47)
        NEW(J,44)=OLD(J,48)
        NEW(J,45)=OLD(J,49)
        NEW(J,46)=OLD(J,50)
        NEW(J,47)=OLD(J,51)
        NEW(J,48)=OLD(J,52)
        NEW(J,49)=OLD(J,5)
        NEW(J,50)=OLD(J,53)
        NEW(J,51)=OLD(J,54)
        NEW(J,52)=OLD(J,6)
        NEW(J,53)=OLD(J,55)
        NEW(J,54)=OLD(J,56)
        NEW(J,55)=OLD(J,57)
        NEW(J,56)=OLD(J,58)
        NEW(J,57)=OLD(J,59)
        NEW(J,58)=OLD(J,60)
        NEW(J,59)=OLD(J,61)
        NEW(J,60)=OLD(J,62)
        NEW(J,61)=OLD(J,7)
        NEW(J,62)=OLD(J,63)
        NEW(J,63)=OLD(J,64)
        NEW(J,64)=OLD(J,8)        

        ELSE
           STOP
        END IF

      ELSE IF (BASE_INFO%TRI_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==3) THEN
        !2D TET LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        ELSE IF(BASE_INFO%B(I)%n==6) THEN
        !2D TET QUAD
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,6)
        NEW(J,6)=OLD(J,5)
        ELSE IF(BASE_INFO%B(I)%n==10) THEN
        !2D TET CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,8)
        NEW(J,7)=OLD(J,9)
        NEW(J,8)=OLD(J,7)
        NEW(J,9)=OLD(J,6)
        NEW(J,10)=OLD(J,10)
        ELSE
          STOP
        END IF

      ELSE IF (BASE_INFO%TET_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==4) THEN
        !3D TET LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==10) THEN
        !3D TET QUAD
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,6)
        NEW(J,7)=OLD(J,7)
        NEW(J,8)=OLD(J,8)
        NEW(J,9)=OLD(J,10)
        NEW(J,10)=OLD(J,9)

        ELSE IF(BASE_INFO%B(I)%n==20) THEN
        !3D TET CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,6)
        NEW(J,7)=OLD(J,7)
        NEW(J,8)=OLD(J,8)
        NEW(J,9)=OLD(J,9)
        NEW(J,10)=OLD(J,10)
        NEW(J,11)=OLD(J,11)
        NEW(J,12)=OLD(J,12)
        NEW(J,13)=OLD(J,15)
        NEW(J,14)=OLD(J,16)
        NEW(J,15)=OLD(J,13)
        NEW(J,16)=OLD(J,14)
        NEW(J,17)=OLD(J,17)
        NEW(J,18)=OLD(J,18)
        NEW(J,19)=OLD(J,19)
        NEW(J,20)=OLD(J,20)
        ELSE
          STOP
        END IF
        
       ELSE
        STOP
       END IF

    END DO



	

END SUBROUTINE ORDER_NUMBERING

! ----------------------------------------------------------------------------------

SUBROUTINE MAKE_UNIQUE
      IMPLICIT NONE

  INTEGER:: I,J,K,L,M,N

! NOW, THE NODE NUMBERING NEEDS TO BE CHANGED FOR ALL QUADRATIC AND CUBIC ELEMENTS


! NOW, SAME NODES NEED MAKE_UNIQUE/IDENTICAL DEFINITION

! M is considered as reference and check V

    MESH_INFO(1)%T=MESH_INFO(1)%T
    MESH_INFO(2)%T=MESH_INFO(2)%T+ArrayOfNodesDefined(1)
    MESH_INFO(3)%T=MESH_INFO(3)%T+ArrayOfNodesDefined(1)+ArrayOfNodesDefined(2)

    IF(ArrayOfNodesDefined(1)==ArrayOfNodesDefined(2)) THEN
    ! copy all node numbers from 2 -> 1
    MESH_INFO(2)%T(:,:)= MESH_INFO(1)%T(:,:)
    ELSE
      IF (BASE_INFO%TRI_BASIS== 1) THEN
        MESH_INFO(2)%T(:,1:3)=MESH_INFO(1)%T(:,1:3)
      ELSE IF (BASE_INFO%TET_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%QUAD_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%HEX_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:8)=MESH_INFO(1)%T(:,1:8)
      ELSE
        STOP
      END IF
    END IF

    IF(ArrayOfNodesDefined(1)==ArrayOfNodesDefined(3)) THEN
    ! copy all node numbers from 2 -> 1
    MESH_INFO(3)%T(:,:)=MESH_INFO(1)%T(:,:)
    ELSE IF(ArrayOfNodesDefined(2)==ArrayOfNodesDefined(3)) THEN
    MESH_INFO(3)%T(:,:)=MESH_INFO(2)%T(:,:)
    ELSE
      IF (BASE_INFO%TRI_BASIS== 1) THEN
        MESH_INFO(3)%T(:,1:3)=MESH_INFO(1)%T(:,1:3)
      ELSE IF (BASE_INFO%TET_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%QUAD_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%HEX_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:8)=MESH_INFO(1)%T(:,1:8)
      ELSE
        STOP
      END IF
    END IF

END SUBROUTINE MAKE_UNIQUE

! ----------------------------------------------------------------------------------

SUBROUTINE READ_NODES

      IMPLICIT NONE

      INTEGER:: I,J,K,L,M,N
      INTEGER:: a,b
      DOUBLE PRECISION,DIMENSION(832,3):: sebo_test_array

      sebo_test_array=0.0

      READ(42,*) NAMz

      OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
      READ(1,*) ArrayOfNodesDefined(1:3)

      TotalNumberOfNodes=ArrayOfNodesDefined(1)+ArrayOfNodesDefined(2)+ArrayOfNodesDefined(3)


! ALLOCATE AND READ MESH NODE INFORMATION

  WRITE(*,*)'Reading Nodes...'

  DO I=1,3
      MESH_INFO(I)%Lx=ArrayOfNodesDefined(I)
      ALLOCATE(MESH_INFO(I)%X(MESH_INFO(I)%Lx,3),STAT=ALLOC_ERROR)
        DO J = 1,MESH_INFO(I)%Lx
          READ(1,*,END=35) MESH_INFO(I)%X(J,1:3)
!	  WRITE(*,*) MESH_INFO(I)%X(J,1:3)
!          READ(1,*,END=35) sebo_test_array(J,1:3)
!	  sebo_test_array(J,1:3)=(/1,2,3/)
        END DO
  END DO




      CLOSE(1)

    ALLOCATE(OPENCMISS_NODE_COORD(TotalNumberOfNodes,3),STAT=ALLOC_ERROR)
    a=1
    b=0


    DO I=1,3
    a=b+1
    b=b+ArrayOfNodesDefined(I)
    OPENCMISS_NODE_COORD(a:b,1:3)=MESH_INFO(I)%X(1:ArrayOfNodesDefined(I),1:3)
    END DO





        RETURN

 35     PRINT *, 'FAILS'
        STOP

  IF(ALLOC_ERROR.NE.0) THEN
    STOP 'Error during allocation'
  END IF

END SUBROUTINE READ_NODES

! ----------------------------------------------------------------------------------

SUBROUTINE READ_ELEMENTS

      IMPLICIT NONE

      INTEGER:: I,J,K,L,M,N

      READ(42,*) NAMz
      CLOSE(42)

      OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
      READ(1,*) NumberOfElementsDefined(1:3)

! ALLOCATE AND READ MESH ELEMENT INFORMATION

      TotalNumberOfNodes=ArrayOfNodesDefined(1)+ArrayOfNodesDefined(2)+ArrayOfNodesDefined(3)


  WRITE(*,*)'Reading Elements...'
  WRITE(*,*)



  DO I=1,3
      MESH_INFO(I)%Lt=NumberOfElementsDefined(I)
      ALLOCATE(MESH_INFO(I)%T(MESH_INFO(I)%Lt,NumberOfNodesPerElement(I)),STAT=ALLOC_ERROR)
        DO J = 1,MESH_INFO(I)%Lt
          READ(1,*,END=30) MESH_INFO(I)%T(J,1:NumberOfNodesPerElement(I))
        END DO
  END DO
      CLOSE(1)

        RETURN

 30     PRINT *, 'FAILS'
        STOP

  IF(ALLOC_ERROR.NE.0) THEN
    STOP 'Error during allocation'
  END IF

END SUBROUTINE READ_ELEMENTS

! ----------------------------------------------------------------------------------

SUBROUTINE PRINT_ON_SCREEN

      IMPLICIT NONE

      INTEGER:: I,J,K,L,M,N

  ! where are the node coordinates stored -> 1 MATRIX



     DO I = 1,TotalNumberOfNodes
               WRITE(*,'("Node ",(I0,4x),1000( F5.3,2x ))')I,OPENCMISS_NODE_COORD(I,1:3)
     END DO


  ! where are the element nodes stored -> 3 MATRICES


    WRITE(*,*)
    WRITE(*,*)
  
    DO I = 1,NumberOfElementsDefined(1)
    WRITE(*,'("M-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_M(I,1:NumberOfNodesPerElement(1))
    END DO

    WRITE(*,*)

    DO I = 1,NumberOfElementsDefined(2)
    WRITE(*,'("V-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_V(I,1:NumberOfNodesPerElement(2))
    END DO

    WRITE(*,*)

    DO I = 1,NumberOfElementsDefined(3)
    WRITE(*,'("P-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_P(I,1:NumberOfNodesPerElement(3))
    END DO


    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)


END SUBROUTINE PRINT_ON_SCREEN

! ----------------------------------------------------------------------------------

END MODULE FLUID_MECHANICS_IO_ROUTINES
