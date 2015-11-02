!> \file
!> \author Ishani Roy & Sander Land
!> \brief This module handles all Monodomain equations routines which use Strang Splitting with hard-coded cell models.
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
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
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

!>This module handles all Monodomain equations routines.
MODULE MONODOMAIN_EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE ELECTROPHYSIOLOGY_CELL_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FITTING_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC MONODOMAIN_CONTROL_LOOP_POST_LOOP
 
  PUBLIC MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP
  
  PUBLIC Monodomain_FiniteElementCalculate

  PUBLIC Monodomain_EquationsSetSolutionMethodSet
  
  PUBLIC Monodomain_EquationsSetSpecificationSet

  PUBLIC Monodomain_ProblemSpecificationSet
  
  PUBLIC MONODOMAIN_EQUATION_PROBLEM_SETUP
  
  PUBLIC MONODOMAIN_PRE_SOLVE,MONODOMAIN_POST_SOLVE
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE MONODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP,TIME_LOOP_PARENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: PARENT_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER,CURRENT_LOOP_ITERATION

    ENTERS("MONODOMAIN_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !Export the dependent field for this time step
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            PROBLEM=>CONTROL_LOOP%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              !Get the solver. 
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)            
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              !Loop over the equations sets associated with the solver
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      NULLIFY(DEPENDENT_REGION)
                      CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                      NULLIFY(PARENT_LOOP)
                      PARENT_LOOP=>CONTROL_LOOP%PARENT_LOOP
                      IF(ASSOCIATED(PARENT_LOOP)) THEN
                        !add the iteration number of the parent loop to the filename
                        NULLIFY(TIME_LOOP_PARENT)
                        TIME_LOOP_PARENT=>PARENT_LOOP%TIME_LOOP
                        IF(ASSOCIATED(TIME_LOOP_PARENT)) THEN
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP_PARENT%OUTPUT_NUMBER
                          CURRENT_LOOP_ITERATION=TIME_LOOP_PARENT%GLOBAL_ITERATION_NUMBER
                          FILENAME="Time_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP_PARENT%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%ITERATION_NUMBER,"*",ERR,ERROR))
                        ELSE
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP%OUTPUT_NUMBER
                          CURRENT_LOOP_ITERATION=TIME_LOOP%GLOBAL_ITERATION_NUMBER
                          FILENAME="Time_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                        ENDIF
                      ELSE
                        OUTPUT_ITERATION_NUMBER=TIME_LOOP%OUTPUT_NUMBER
                        CURRENT_LOOP_ITERATION=TIME_LOOP%GLOBAL_ITERATION_NUMBER
                        FILENAME="Time_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                          & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                      ENDIF
                      METHOD="FORTRAN"
                      IF(OUTPUT_ITERATION_NUMBER/=0.AND.MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0) THEN
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index "// &
                        & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                        & " in the solver mapping."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ELSE
                  CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE DEFAULT
          LOCAL_ERROR="The control loop type of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MONODOMAIN_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MONODOMAIN_CONTROL_LOOP_POST_LOOP
  !
  !================================================================================================================================
  !
  !

  !>Sets the problem specification for a monodomain equation set class.
  SUBROUTINE Monodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Monodomain_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<3) THEN
        CALL FlagError("Equations set specification must have at least three entries for a monodomain class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        subtype=specification(3)
        SELECT CASE(subtype)
        CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVstring(subtype,"*",err,error))// &
            & " is not valid for a Monodomain equation type of a Strang splitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Set full specification
        IF(ALLOCATED(equationsSet%specification)) THEN
          CALL FlagError("Equations set specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(equationsSet%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
        END IF
        equationsSet%specification(1:3)=[EQUATIONS_SET_BIOELECTRICS_CLASS, &
          & EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE,subtype]
      CASE DEFAULT
        localError="Equations set equation type "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for a monodomain equations set class."
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

    CALL Exits("Monodomain_EquationsSetSpecificationSet")
    RETURN
999 CALL Errors("Monodomain_EquationsSetSpecificationSet",err,error)
    CALL Exits("Monodomain_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a monodomain problem class.
  SUBROUTINE Monodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType,problemSubtype

    CALL Enters("Monodomain_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=3) THEN
        problemType=problemSpecification(2)
        SELECT CASE(problemType)
        CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
          problemSubtype=problemSpecification(3)
          SELECT CASE(problemSubtype)
          CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
            & PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
            !ok
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
              & " is not valid for a Monodomain equation type of a Strang splitting  problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          IF(ALLOCATED(problem%specification)) THEN
            CALL FlagError("Problem specification is already allocated.",err,error,*999)
          ELSE
            ALLOCATE(problem%specification(3),stat=err)
            IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
          END IF
          problem%specification(1:3)=[PROBLEM_BIOELECTRICS_CLASS,PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE,problemSubtype]
        CASE DEFAULT
          localError="Problem equation type "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a monodomain problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Monodomain problem specification must have a type.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated",err,error,*999)
    END IF

    CALL Exits("Monodomain_ProblemSpecificationSet")
    RETURN
999 CALL Errors("Monodomain_ProblemSpecificationSet",err,error)
    CALL Exits("Monodomain_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Monodomain equation finite element equations set.
  SUBROUTINE Monodomain_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nh,nhs,ni,ns,nj
    REAL(DP) :: RWG,SUM,Df, Dt, D(3,3), f(3), fnorm
    REAL(DP) :: DPHIDX(3,8) ! assumes <= 8 basis functions / DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS <= 8
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: DAMPING_MATRIX,STIFFNESS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    ENTERS("Monodomain_FiniteElementCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
          CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
            & err,error,*999)
        END IF
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE,EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)

          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR

          IF(.NOT.(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX.OR.RHS_VECTOR%UPDATE_VECTOR)) RETURN ! no updates -> return immediately
 
          ! set up opencmiss interpolation and basis
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)


          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            ! get interpolated geometric and material interpolated point
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &

              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

            ! calculate weight = det J * gauss pt weight
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

            ! basis function chain rule taken out of the inner loop.
            DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
              DPHIDX(nj,ms)=0.0_DP
              DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                DPHIDX(nj,ms)=DPHIDX(nj,ms) + &
                             & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                             & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(ni,nj)
              ENDDO !ni
            ENDDO !nj
            ENDDO !ms

            ! Diffusion tensor  D  =  Dt I + (Df - Dt) f f^T  where Dt and Df are diffusivity/conductivity in fiber/transverse directions
            Df = EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV) ! 2 = Df
            Dt = EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV) ! 3 = Dt
            fnorm = 0.0
            DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
              f(nj) = EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3+nj,NO_PART_DERIV) ! 4,5[,6] = f
              fnorm = fnorm + f(nj)*f(nj)
            ENDDO
            ! normalize f, and fill in default for 0,0,0 -> 1,0,0
            fnorm = SQRT(fnorm)
            IF(fnorm < 1e-6) THEN
              f = (/ 1.0, 0.0, 0.0 /) ! default
            ELSE
              f = f / fnorm
            ENDIF
            DO ni=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
              D(ni,:)  = 0.0
              D(ni,ni) = Dt
              DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                D(ni,nj) = D(ni,nj) + (Df - Dt) * f(ni) * f(nj)
              ENDDO
            ENDDO

            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns. TODO: use symmetry?
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        SUM=0.0_DP
                        DO ni=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS 
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          SUM=SUM + D(ni,nj) * DPHIDX(ni,ms) * DPHIDX(nj,ns)
                        ENDDO !nj
                        ENDDO !ni
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG ! Aij = int D_ij * dphi_m/dx_i * dphi_n/dx_j 
                      ENDIF
 
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN ! non mass lumped version
                         DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG ! int phi_m phi_n
                      ENDIF 

                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
!                IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN  ! mass lumnped version
!                  DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,mhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,mhs)+ &
!                  & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*RWG   !  // int phi_m
!                ENDIF
              ENDDO !ms
            ENDDO !mh
          IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP 
          ENDDO !ng
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN ! non mass lumped version
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                           & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                           & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                      ENDIF 
                    ENDDO !ns
                  ENDDO !nh
                ENDIF

                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh) 
 !               IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN ! mass lumped version
 !                  DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
 !                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
 !               ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Monodomain equation type of a Strang splitting equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("Monodomain_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Monodomain_FiniteElementCalculate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Monodomain_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the Monodomain equation type of a Strang splitting equations set class.
  SUBROUTINE MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Monodomain equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE,EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        CALL Monodomain_EquationsSetSubtypeSetUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)        
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation type of a Strang splitting equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Monodomain equation type of an Strang splitting equations set class.
  SUBROUTINE Monodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE,EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation type of monodomain Strang splitting equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Monodomain_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Monodomain_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("Monodomain_EquationsSetSolutionMethodSet")
    RETURN 1

  END SUBROUTINE Monodomain_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the Monodomain equation.
  SUBROUTINE Monodomain_EquationsSetSubtypeSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS, NUM_COMP
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_EQUATION_SET_SUBTYPE_SETUP",ERR,ERROR,*999)
 
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)


    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE .OR. &
          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Monodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999) 

          

              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,&
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,&
                & 1 ,ERR,ERROR,*999)

              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid or not implemented"
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              ! user specified field
              CALL FlagError("No user specified field supported!",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation"
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
              CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
                 NUM_COMP = 4
              CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
                 NUM_COMP = 19
              CASE DEFAULT
                CALL FlagError("Invalid cell model equations set subtype",ERR,ERROR,*999)
              END SELECT    

              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                & INDEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,&
                   & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,&
                   &(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
         
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
                   &NUM_COMP,ERR,ERROR,*999)
              CALL FIELD_DOF_ORDER_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
              &  FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,ERR,ERROR,*999) ! dofs continuous, so first + (x-1) is x'th component index

              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)

              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              END SELECT
            ELSE
              ! user specified field
              CALL FlagError("No user specified field supported!",ERR,ERROR,*999)
            ENDIF

          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
              CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
                CALL BUENO_OROVIO_INITIALIZE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999) ! initialize to y0
              CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
                CALL TENTUSSCHER06_INITIALIZE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999) ! initialize to y0 
              CASE DEFAULT
                CALL FlagError("Invalid cell model equations set subtype",ERR,ERROR,*999)
              END SELECT              
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation"
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
           CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
             IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
               IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  !Materials field components are
                  ! 1. activation  factor (usually 0.0 or 1.0)
                  ! 2,3 for fiber/transverse conductivity   . defaults to constant interpolation 
                  ! 4,5[,6] : fiber unit vector in dimension
                  ! 7: out - activation times
                NUMBER_OF_MATERIALS_COMPONENTS= 7 !NUMBER_OF_DIMENSIONS + 3
                 !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)

                ! 1st = activation = node based, 2 3 diffusion constants
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 3,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ! 4 5 (6) fiber unit vector
                DO component_idx=1,3 !NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx+3,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx+3,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx

                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)

                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                ! user specified field
                CALL FlagError("No user specified field supported!",ERR,ERROR,*999)
              ENDIF
             ELSE

               CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
             ENDIF
           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !Materials field components are 1 for each dimension 
                !i.e., k in div(k.grad(u(x)))
                NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS                             
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
! ! ! Upto here
           CASE DEFAULT
             LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
               & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
               & " is invalid for a Monodomain equation."
             CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
           END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              !!! Check this 
              !CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
              !CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
               ! & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              !!!! Check the above two lines
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                !CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                 ! & ERR,ERROR,*999)
                !CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                 ! & ERR,ERROR,*999)
                 CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & (/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                    (/EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE/),ERR,ERROR,*999)     
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
            CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid or not implemented."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a Monodomain equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a Monodomain equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Monodomain_EquationsSetSubtypeSetup")
    RETURN
999 ERRORSEXITS("Monodomain_EquationsSetSubtypeSetup",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSubtypeSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the Monodomain solution.
  SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Monodomain equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a monodomain equation problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        CALL Monodomain_ProblemStrangSplittingSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation Strang splitting problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SETUP
  
!
  !================================================================================================================================
  !
 
  !>Sets up the Monodomain solution.
  SUBROUTINE Monodomain_ProblemStrangSplittingSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Monodomain equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Monodomain_ProblemStrangSplittingSetup",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a monodomain Strang splitting problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE,PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        CALL MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation type of a Strang splitting problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Monodomain_ProblemStrangSplittingSetup")
    RETURN
999 ERRORS("Monodomain_ProblemStrangSplittingSetup",ERR,ERROR)
    EXITS("Monodomain_ProblemStrangSplittingSetup")
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemStrangSplittingSetup
  
  !
  !================================================================================================================================
  !

  !>Sets up the linear Monodomain equations solution.
  SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a monodomain equation problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE.OR. &
          & PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
           !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a first order dynamic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
           SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a Monodomain equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a Monodomain equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a monodomain problem class.
  SUBROUTINE MONODOMAIN_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD, INDEPENDENT_FIELD

    ENTERS("MONODOMAIN_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must at least two entries for a monodomain equation problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        DEPENDENT_FIELD=>SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
          & EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
        INDEPENDENT_FIELD=>SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
          & EQUATIONS_SETS(1)%PTR%INDEPENDENT%INDEPENDENT_FIELD

        CALL Field_ParametersToFieldParametersCopy(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 1,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,ERR,ERROR,*999)
        CALL Field_ParametersToFieldParametersCopy(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 1,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,ERR,ERROR,*999) ! also to prev.

      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a monodomain problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("MONODOMAIN_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a monodomain problem class.
  SUBROUTINE MONODOMAIN_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP) :: TMPV, TMPA
    INTEGER(INTG) :: I
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD, MATERIAL_FIELD,  INDEPENDENT_FIELD, GEOMETRIC_FIELD

    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: DAMPING_MATRIX,STIFFNESS_MATRIX
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES

    ENTERS("MONODOMAIN_POST_SOLVE",ERR,ERROR,*999)
  
    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a monodomain equation problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) ! CAN NOT GET EQN TYPE?! 


       ! Don't update in next step
        EQUATIONS=>SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%EQUATIONS
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
        DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        STIFFNESS_MATRIX%UPDATE_MATRIX = .FALSE.
        DAMPING_MATRIX%UPDATE_MATRIX   = .FALSE.
        RHS_VECTOR%UPDATE_VECTOR       = .FALSE.

       ! integrate   cell models
        GEOMETRIC_FIELD => SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                           & EQUATIONS_SETS(1)%PTR%GEOMETRY%GEOMETRIC_FIELD
        DEPENDENT_FIELD => SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                           & EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
        MATERIAL_FIELD  => SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                           & EQUATIONS_SETS(1)%PTR%MATERIALS%MATERIALS_FIELD
        INDEPENDENT_FIELD => SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                           & EQUATIONS_SETS(1)%PTR%INDEPENDENT%INDEPENDENT_FIELD

        CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & 1, INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1,ERR,ERROR,*999) ! dependent -> independent

        SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE)
          CALL BUENO_OROVIO_INTEGRATE(INDEPENDENT_FIELD,MATERIAL_FIELD,&
          &  CONTROL_LOOP%TIME_LOOP%CURRENT_TIME-CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT, CONTROL_LOOP%TIME_LOOP%CURRENT_TIME&  ! from t-dt to t
          & ,ERR,ERROR,*999)
        CASE(PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
           CALL TENTUSSCHER06_INTEGRATE(INDEPENDENT_FIELD,MATERIAL_FIELD,&
          &  CONTROL_LOOP%TIME_LOOP%CURRENT_TIME-CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT, CONTROL_LOOP%TIME_LOOP%CURRENT_TIME&  ! from t-dt to t
          & ,ERR,ERROR,*999)
        CASE DEFAULT
          CALL FlagError("Invalid cell model subtype",ERR,ERROR,*999)
        END SELECT

        DO I=1,INDEPENDENT_FIELD%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY%NODES%NUMBER_OF_NODES
          !Default to version 1 of each derivative
          CALL FIELD_PARAMETER_SET_GET_NODE(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,I,1,TMPV,&  ! get local node?
               & ERR,ERROR,*999)
          IF (TMPV > 0) THEN
           !Default to version 1 of each derivative
           CALL FIELD_PARAMETER_SET_GET_NODE(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,I,7,TMPA,&
                & ERR,ERROR,*999) 
            IF(ABS(TMPA)<ZERO_TOLERANCE) THEN
              !Default to version 1 of each derivative
              CALL FIELD_PARAMETER_SET_UPDATE_NODE(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,I,7,&
                   &CONTROL_LOOP%TIME_LOOP%CURRENT_TIME, ERR,ERROR,*999)  
            ENDIF      
          ENDIF
        ENDDO

        IF(MOD(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME+1e-6,5.0)<1e-3) THEN
          WRITE(*,*), 'T=',CONTROL_LOOP%TIME_LOOP%CURRENT_TIME
        ENDIF

      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a monodomain problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_POST_SOLVE")
    RETURN
999 ERRORSEXITS("MONODOMAIN_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_POST_SOLVE

  !
  !================================================================================================================================
  !

END MODULE MONODOMAIN_EQUATIONS_ROUTINES
