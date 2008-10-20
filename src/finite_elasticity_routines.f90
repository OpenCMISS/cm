!> \file
!> $Id: finite_elasticity_routines.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief This module handles all finite elasticity routines.
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
!> Contributor(s): Kumar Mithraratne
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

!>This module handles all finite elasticity routines.
MODULE FINITE_ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS  
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE,FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SETUP,FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET, &
    & FINITE_ELASTICITY_PROBLEM_SETUP
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a finite elasticity finite element equations set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
     
    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual and RHS vectors for a finite elasticity finite element equations set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: field_interpolation_parameters_undef, & 
      & field_interpolation_parameters_def,field_interpolation_parameters_fibre, &
      & field_interpolation_parameters_material,field_interpolation_parameters_hydrostatic        
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolated_point_hydrostatic,interpolated_point_undef, &
      &interpolated_point_def,interpolated_point_fibre,interpolated_point_material 
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,undef_field,def_field,fibre_field,material_field,hydrostatic_field             
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG) :: dof_idx,gauss_point_idx,NUMBER_OF_DOFS,NEXT_NUMBER,NXI,xi_idx,NGT,nj_idx,ns,i,j,TOTAL_DOFS
    INTEGER(INTG), ALLOCATABLE :: NGXI(:)  
    REAL(DP) :: CAUCHY_TENSOR(3,3),DZDNU(3,3),Jznu,Jxxi,WG
    REAL(DP), ALLOCATABLE :: DFDZ(:,:),RE(:),XI(:),XIG(:,:)  
       
    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)  
	!WRITE(*,*) ' Calculating Residulas *****'  
	
	REGION=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%REGION
	DECOMPOSITION=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%DECOMPOSITION	
	GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
	
	!undef field	
        undef_field=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(undef_field,1,field_interpolation_parameters_undef,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_INITIALISE(field_interpolation_parameters_undef,interpolated_point_undef,ERR,ERROR,*999)
	CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(1,ELEMENT_NUMBER,field_interpolation_parameters_undef,ERR,ERROR,*999)        
        !QUADRATURE_SCHEME_UNDEF=field_interpolation_parameters_undef%BASES(1)%PTR%QUADRATURE%TYPE      	
	NXI=field_interpolation_parameters_undef%BASES(1)%PTR%NUMBER_OF_XI
	ALLOCATE(NGXI(NXI))
	NGT=1
	DO xi_idx=1,NXI,1
	  NGXI(xi_idx)=field_interpolation_parameters_undef%BASES(1)%PTR%QUADRATURE%NUMBER_OF_GAUSS_XI(xi_idx)
	  NGT=NGT*NGXI(xi_idx)
	ENDDO
	ALLOCATE(XIG(NXI,NGT))
	ALLOCATE(XI(NXI))
 
        !def field
	CALL FIELD_NEXT_NUMBER_FIND(REGION,NEXT_NUMBER,ERR,ERROR,*999)
        CALL FIELD_CREATE_START(NEXT_NUMBER,REGION,def_field,ERR,ERROR,*999)
	def_field%DECOMPOSITION=>DECOMPOSITION
	def_field%GEOMETRIC_FIELD=>GEOMETRIC_FIELD
        CALL FIELD_CREATE_FINISH(undef_field%REGION,def_field,ERR,ERROR,*999)
	NUMBER_OF_DOFS=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(1)%NUMBER_OF_DOFS			
        DO dof_idx=1,NUMBER_OF_DOFS,1
	  def_field%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS%CMISS%DATA_DP(dof_idx)= &
	    & undef_field%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS%CMISS%DATA_DP(dof_idx)+ &
	    & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%PARAMETER_SETS%PARAMETER_SETS(2)%PTR%PARAMETERS%CMISS%DATA_DP(dof_idx) 
	                                             !??? should be %PARAMETER_SETS(1) - kmith
	ENDDO	
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(def_field,1,field_interpolation_parameters_def,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_INITIALISE(field_interpolation_parameters_def,interpolated_point_def,ERR,ERROR,*999)
	CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(1,ELEMENT_NUMBER,field_interpolation_parameters_def,ERR,ERROR,*999)        
         
	!hydrostatic pressure field
	CALL FIELD_NEXT_NUMBER_FIND(REGION,NEXT_NUMBER,ERR,ERROR,*999)                     !get next field number
	CALL FIELD_CREATE_START(NEXT_NUMBER,REGION,hydrostatic_field,ERR,ERROR,*999)       !create field called hydrostatic_field
	CALL FIELD_TYPE_SET(NEXT_NUMBER,REGION,FIELD_GENERAL_TYPE,ERR,ERROR,*999)          !set type to FIELD_GENERAL_TYPE : see field_routine constants
	CALL FIELD_MESH_DECOMPOSITION_SET(NEXT_NUMBER,REGION,DECOMPOSITION,ERR,ERROR,*999) !associate a decomposition (geometric_field decomposition)
	CALL FIELD_GEOMETRIC_FIELD_SET(NEXT_NUMBER,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)  !associate a geomtery (geomteric_field)
	CALL FIELD_NUMBER_OF_COMPONENTS_SET(NEXT_NUMBER,REGION,1,ERR,ERROR,*999)           !set no. of field components to 1, -i.e. scalar 
        CALL FIELD_COMPONENT_INTERPOLATION_SET(NEXT_NUMBER,1,1,REGION, &                   !set interpolation to FIELD_ELEMENT_BASED_INTERPOLATION : 
	  & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)                              !see field_routine constants 	
	CALL FIELD_CREATE_FINISH(undef_field%REGION,hydrostatic_field,ERR,ERROR,*999)
	hydrostatic_field%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS%CMISS%DATA_DP(1)=-5.0_DP !set the value to -5.0 : initial guess.
	  
        !fibre field
	fibre_field=>EQUATIONS_SET%GEOMETRY%FIBRE_FIELD
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(fibre_field,1,field_interpolation_parameters_fibre,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_INITIALISE(field_interpolation_parameters_fibre,interpolated_point_fibre,ERR,ERROR,*999)
	CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(1,ELEMENT_NUMBER,field_interpolation_parameters_fibre,ERR,ERROR,*999)        
		
        !material field	
	material_field=>EQUATIONS_SET%MATERIALS%MATERIAL_FIELD
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(material_field,1,field_interpolation_parameters_material,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINT_INITIALISE(field_interpolation_parameters_material,interpolated_point_material,ERR,ERROR,*999)
	CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(1,ELEMENT_NUMBER,field_interpolation_parameters_material,ERR,ERROR,*999)        
        
	TOTAL_DOFS=undef_field%VARIABLES(1)%NUMBER_OF_DOFS+hydrostatic_field%VARIABLES(1)%NUMBER_OF_DOFS
	ALLOCATE(RE(TOTAL_DOFS)) 
	ns=0              	    
        DO nj_idx=1,3,1
	  DO dof_idx=1,undef_field%VARIABLES(1)%COMPONENTS(nj_idx)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS,1
	    ns=ns+1
	    RE(ns)=0.0_DP
	    !EQUATIONS_SET%EQUATIONS%EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR(ns)=0.0_DP
	    DO gauss_point_idx=1,NGT,1          
	      !WRITE(*,'(2(i5))') ELEMENT_NUMBER,gauss_point_idx  		  
	      DO xi_idx=1,NXI,1		
	        XIG(xi_idx,gauss_point_idx)=interpolated_point_undef%INTERPOLATION_PARAMETERS%BASES(1)%PTR%QUADRATURE%SCHEMES(1)% &
	          & PTR%GAUSS_POSITIONS(xi_idx,gauss_point_idx)	        	
	        XI(xi_idx)=XIG(xi_idx,gauss_point_idx)	    	  
	      ENDDO	
	      WG=interpolated_point_undef%INTERPOLATION_PARAMETERS%BASES(1)%PTR%QUADRATURE%SCHEMES(1)%PTR%GAUSS_WEIGHTS(gauss_point_idx)
	      
	      CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_undef,ERR,ERROR,*999)
	      CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_def,ERR,ERROR,*999)	  
	      CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_fibre,ERR,ERROR,*999)	  
	      CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_material,ERR,ERROR,*999)
		
	      CALL FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR(interpolated_point_def,interpolated_point_fibre, &
	        & interpolated_point_undef,ELEMENT_NUMBER,gauss_point_idx,DZDNU,Jxxi,Jznu,ERR,ERROR,*999)    
	
              CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(hydrostatic_field,interpolated_point_material,CAUCHY_TENSOR,DZDNU,ERR,ERROR,*999)
		         	  
	      CALL FINITE_ELASTICITY_GAUSS_DFDZ(interpolated_point_def,DFDZ,XI,ERR,ERROR,*999)	        	  
	  	  	  
	      i=nj_idx
	      DO j=1,3,1
		!EQUATIONS_SET%EQUATIONS%EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR(ns)= &
		!  & EQUATIONS_SET%EQUATIONS%EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR(ns)+ &
		!  & WG*Jxxi*Jznu*(CAUCHY_TENSOR(i,j)*DFDZ(dof_idx,j))	
		RE(ns)=RE(ns)+WG*Jxxi*Jznu*(CAUCHY_TENSOR(i,j)*DFDZ(dof_idx,j))   
	      ENDDO				
	  	  	      
	      DEALLOCATE(DFDZ)    	             
            ENDDO !gauss_idx
		
	  ENDDO !dof_idx	
	ENDDO !nj_idx
	
	ns=19
	WRITE(*,'(i5,f12.7)') ns,RE(ns) !EQUATIONS_SET%EQUATIONS%EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR(ns)
		
	!Hydrostatic pressure equation/s
	DO dof_idx=ns+1,TOTAL_DOFS,1	    
	  ns=ns+1
	  RE(ns)=0.0_DP
	   
	  DO gauss_point_idx=1,NGT,1          	  
	    DO xi_idx=1,NXI,1		
	      XIG(xi_idx,gauss_point_idx)=interpolated_point_undef%INTERPOLATION_PARAMETERS%BASES(1)%PTR%QUADRATURE%SCHEMES(1)% &
	        & PTR%GAUSS_POSITIONS(xi_idx,gauss_point_idx)	        	
	      XI(xi_idx)=XIG(xi_idx,gauss_point_idx)	    	  
	    ENDDO	
	    WG=interpolated_point_undef%INTERPOLATION_PARAMETERS%BASES(1)%PTR%QUADRATURE%SCHEMES(1)%PTR%GAUSS_WEIGHTS(gauss_point_idx)
	      
	    CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_undef,ERR,ERROR,*999)
	    CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_def,ERR,ERROR,*999)	  
	    CALL FIELD_INTERPOLATE_XI(PART_DERIV_S1,XI,interpolated_point_fibre,ERR,ERROR,*999)	  
		
	    CALL FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR(interpolated_point_def,interpolated_point_fibre, &
	      & interpolated_point_undef,ELEMENT_NUMBER,gauss_point_idx,DZDNU,Jxxi,Jznu,ERR,ERROR,*999)    
	  
	    RE(ns)=RE(ns)+WG*Jxxi*1.0_DP*(Jznu-1.0_DP)
	    
	  ENDDO
	ENDDO
	!WRITE(*,'(i5,f14.9)') ns,RE(25)
	
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !
  
  !>Evaluates the deformation gradient tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR(INTERPOLATED_POINT_DEF,INTERPOLATED_POINT_FIBRE, &
    & INTERPOLATED_POINT_UNDEF,ELEMENT_NUMBER,GAUSS_PT_NUMBER,DZDNU,Jxxi,Jznu,ERR,ERROR,*)    

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT_DEF,INTERPOLATED_POINT_FIBRE,INTERPOLATED_POINT_UNDEF    
    INTEGER(INTG) :: ELEMENT_NUMBER,GAUSS_PT_NUMBER      
    REAL(DP) :: DZDNU(3,3),Jxxi,Jznu    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string   
    !Local Variables
    INTEGER(INTG) :: derivative_idx,nj_idx,xi_idx 
    REAL(DP) :: DNUDX(3,3),DNUDXI(3,3),DXDNU(3,3),DXIDNU(3,3),DXDXI(3,3),DZDXI(3,3),Jnuxi
    INTEGER(INTG) :: i,j    

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR",ERR,ERROR,*999)
    	 
    DO nj_idx=1,3,1
      DO xi_idx=1,3,1
        derivative_idx=INT((xi_idx*xi_idx+xi_idx+2+0.01_DP)/2.0_DP)  !2,4,7      
        DXDXI(nj_idx,xi_idx)=interpolated_point_undef%VALUES(nj_idx,derivative_idx)  !dx/dxi
        DZDXI(nj_idx,xi_idx)=interpolated_point_def%VALUES(nj_idx,derivative_idx) !dz/dxi
      ENDDO	  
    ENDDO 

    CALL FINITE_ELASTICITY_GAUSS_DXDNU(INTERPOLATED_POINT_FIBRE,DXDNU,DXDXI,ERR,ERROR,*999)

    CALL MATRIX_TRANSPOSE(DXDNU,DNUDX,ERR,ERROR,*999) !dx/dnu is orthogonal. Therefore transpose is inverse
          	  
    CALL MATRIX_PRODUCT(DNUDX,DXDXI,DNUDXI,ERR,ERROR,*999) ! dnu/dxi = dnu/dx * dx/dxi
    CALL INVERT(DNUDXI,DXIDNU,Jnuxi,ERR,ERROR,*999) ! dxi/dnu 
	  
    CALL MATRIX_PRODUCT(DZDXI,DXIDNU,DZDNU,ERR,ERROR,*999) ! dz/dnu = dz/dxi * dxi/dnu  (deformation gradient tensor, F)	
            
    Jxxi=DETERMINANT(DXDXI,ERR,ERROR)
    Jznu=DETERMINANT(DZDNU,ERR,ERROR)
              
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADEINT_TENSOR

  !
  !================================================================================================================================
  !
  
  !>Evaluates dx/dnu(undeformed-material cs) tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DXDNU(INTERPOLATED_POINT,DXDNU,DXDXI,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP) :: DXDNU(:,:),DXDXI(:,:)    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: fibre_idx,i,j,nj_idx 
    INTEGER(INTG) :: vector(3) = (/1,2,3/)
    REAL(DP) :: A(3),ANGLE(3),B(3),C(3),DXDNU1(3,3),DXDNU2(3,3),DXDNU3(3,3),R(3,3), &
      & R1(3,3),R2(3,3),R3(3,3)

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_DXDNU",ERR,ERROR,*999)
        
    !First calculate default dx/dnu
    DO nj_idx=1,3,1 	
      A(nj_idx)=DXDXI(nj_idx,1) !default material_1 dir.
    ENDDO 
    CALL CROSS_PRODUCT(DXDXI(vector,1),DXDXI(vector,2),C,ERR,ERROR,*999) !default material_3 dir.    
    CALL CROSS_PRODUCT(C,A,B,ERR,ERROR,*999) ! default material_2 dir.      
    
    DO nj_idx=1,3,1 	
      DXDNU(nj_idx,1)=A(nj_idx)/L2NORM(A) 
      DXDNU(nj_idx,2)=B(nj_idx)/L2NORM(B)      
      DXDNU(nj_idx,3)=C(nj_idx)/L2NORM(C)        
    ENDDO 
    
    !The normalised DXDNU contains the transformation(rotaion) from spatial CS -> material CS 
    DO i=1,3,1
      DO j=1,3,1
        R(i,j)=DXDNU(i,j) 
      ENDDO
    ENDDO
    
    !Now transform(rotate) the default material CS to align with spatial CS.
    !This does not have to be done. Since both CSs are ortonomral, v1-v2-v3 
    !end up having coordinates (1,0,0), (0,1,0) and (0,0,1)
    !Then rotate by angles specified (old CMISS - ipfibr file values)     
    DO fibre_idx=1,3,1     
      ANGLE(fibre_idx)=INTERPOLATED_POINT%VALUES(fibre_idx,1)
    ENDDO
     
    DO i=1,3,1
      DO j=1,3,1
        R1(i,j)=0.0_DP 
        R2(i,j)=0.0_DP 
        R3(i,j)=0.0_DP 
	IF (i==j) THEN
	  R1(i,j)=1.0_DP
	  R2(i,j)=1.0_DP	  
	  R3(i,j)=1.0_DP
	ENDIF  	  		
      ENDDO
    ENDDO
        
    R3(1,1)=cos(ANGLE(1))    !angles are in radians
    R3(1,2)=-sin(ANGLE(1))
    R3(2,1)=sin(ANGLE(1))
    R3(2,2)=cos(ANGLE(1))
    R2(1,1)=cos(ANGLE(2))
    R2(1,3)=sin(ANGLE(2))
    R2(3,1)=-sin(ANGLE(2))
    R2(3,3)=cos(ANGLE(2))
    R1(2,2)=cos(ANGLE(3))
    R1(2,3)=-sin(ANGLE(3))
    R1(3,2)=sin(ANGLE(3))
    R1(3,3)=cos(ANGLE(3))
    
    CALL MATRIX_PRODUCT(R3,DXDNU,DXDNU3,ERR,ERROR,*999)   !rotate about v3 => v1'-v2'-v3
    CALL MATRIX_PRODUCT(R2,DXDNU3,DXDNU2,ERR,ERROR,*999)  !rotate about v2' => v1''-v2'-v3'     
    CALL MATRIX_PRODUCT(R1,DXDNU2,DXDNU1,ERR,ERROR,*999)  !rotate about v1'' => v1''-v2''-v3''
    
    !Inverse-rotate v1''-v2''-v3''
    CALL MATRIX_PRODUCT(R,DXDNU1,DXDNU,ERR,ERROR,*999)  
            
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DXDNU")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_DXDNU",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DXDNU")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_DXDNU

  !
  !================================================================================================================================
  !

  !>Evaluates the Cauchy stress tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(FIELD,INTERPOLATED_POINT,CAUCHY_TENSOR,DZDNU,ERR,ERROR,*)

    !Argument variables    
    TYPE(FIELD_TYPE), POINTER :: FIELD    
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT 
    REAL(DP), INTENT(OUT) :: CAUCHY_TENSOR(:,:)
    REAL(DP), INTENT(IN) ::  DZDNU(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j    
    REAL(DP) :: AZL(3,3),C(0:1,0:1),DZDNUT(3,3),PIOLA_TENSOR(3,3),TEMP(3,3),Jznu,P 
        
    CALL ENTERS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR,*999)
	
    C(1,0)=INTERPOLATED_POINT%VALUES(1,1)
    C(0,1)=INTERPOLATED_POINT%VALUES(2,1)
    P=FIELD%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS%CMISS%DATA_DP(1)
    
    CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)    
    CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999) ! AZL = F'*F (deformed covariant or right cauchy defromation tensor, C)	    
        
    PIOLA_TENSOR(1,1)=2.0_DP*(C(1,0)+C(0,1)*(AZL(2,2)+AZL(3,3))+P*(AZL(2,2)*AZL(3,3)-AZL(2,3)*AZL(3,2)))
    PIOLA_TENSOR(1,2)=2.0_DP*(       C(0,1)*(-AZL(2,1))        +P*(AZL(2,3)*AZL(3,1)-AZL(2,1)*AZL(3,3)))
    PIOLA_TENSOR(1,3)=2.0_DP*(       C(0,1)*(-AZL(3,1))        +P*(AZL(2,1)*AZL(3,2)-AZL(2,2)*AZL(3,1)))            
    PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
    PIOLA_TENSOR(2,2)=2.0_DP*(C(1,0)+C(0,1)*(AZL(3,3)+AZL(1,1))+P*(AZL(1,1)*AZL(3,3)-AZL(1,3)*AZL(3,1)))
    PIOLA_TENSOR(2,3)=2.0_DP*(       C(0,1)*(-AZL(3,2))        +P*(AZL(1,2)*AZL(3,1)-AZL(1,1)*AZL(3,2)))    
    PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)       
    PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3) 
    PIOLA_TENSOR(3,3)=2.0_DP*(C(1,0)+C(0,1)*(AZL(1,1)+AZL(2,2))+P*(AZL(1,1)*AZL(2,2)-AZL(1,2)*AZL(2,1)))    

    Jznu=DETERMINANT(DZDNU,ERR,ERROR)

    CALL MATRIX_PRODUCT(DZDNU,PIOLA_TENSOR,TEMP,ERR,ERROR,*999)     
    CALL MATRIX_PRODUCT(TEMP,DZDNUT,CAUCHY_TENSOR,ERR,ERROR,*999)     

    DO i=1,3,1
      DO j=1,3,1
        CAUCHY_TENSOR(i,j)=CAUCHY_TENSOR(i,j)/Jznu
      ENDDO
    ENDDO
      
    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR


  !
  !================================================================================================================================
  !

  !>Evaluates df/dz (derivative of interpolation function wrt deformed coord) matrix at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ(INTERPOLATED_POINT,DFDZ,XI,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT       
    REAL(DP), ALLOCATABLE :: DFDZ(:,:),XI(:)     
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string    
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG) :: derivative_idx,NJT,nj_idx,nk_idx,NKT,nn_idx,NNT,ns,NST,NXIT,xi_idx 
    REAL(DP), ALLOCATABLE :: DFDXI(:,:),DXIDZ(:,:),DZDXI(:,:)
    REAL(DP) :: Jzxi
 
    CALL ENTERS("FINITE_ELASTICITY_GAUSS_DFDZ",ERR,ERROR,*999)

    BASIS=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%BASES(1)%PTR
    NJT=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
    NNT=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%BASES(1)%PTR%NUMBER_OF_NODES
    NKT=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%BASES(1)%PTR%MAXIMUM_NUMBER_OF_DERIVATIVES
    NXIT=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%BASES(1)%PTR%NUMBER_OF_XI
    !NST=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%BASES(1)%PTR%NUMBER_OF_ELEMENT_PARAMETERS
    
    ALLOCATE(DXIDZ(NXIT,NJT))
    ALLOCATE(DZDXI(NJT,NXIT))    
    NST=0
    DO nn_idx=1,NNT,1
      DO nk_idx=1,NKT,1
        NST=NST+1    
      ENDDO
    ENDDO
    ALLOCATE(DFDZ(NST,NJT))
    ALLOCATE(DFDXI(NST,NXIT))

    ns=0
    DO nn_idx=1,NNT,1
      DO nk_idx=1,NKT,1
        ns=ns+1    
	DO xi_idx=1,NXIT,1
          derivative_idx=INT((xi_idx*xi_idx+xi_idx+2+0.01_DP)/2.0_DP)  !2,4,7 
          DFDXI(ns,xi_idx)=BASIS_LHTP_BASIS_EVALUATE(BASIS,nn_idx,nk_idx,derivative_idx,XI,ERR,ERROR) !df/dxi 
	ENDDO
      ENDDO
    ENDDO

    DO nj_idx=1,NJT,1
      DO xi_idx=1,NXIT,1
        derivative_idx=INT((xi_idx*xi_idx+xi_idx+2+0.01_DP)/2.0_DP)  !2,4,7      
        DZDXI(nj_idx,xi_idx)=INTERPOLATED_POINT%VALUES(nj_idx,derivative_idx) !dz/dxi
      ENDDO	  
    ENDDO 
    CALL INVERT(DZDXI,DXIDZ,Jzxi,ERR,ERROR,*999) !dxi/dz 
    
    DO nj_idx=1,NJT,1
      ns=0
      DO nn_idx=1,NNT,1
        DO nk_idx=1,NKT,1
          ns=ns+1 
    	  DFDZ(ns,nj_idx)=0.0_DP 
    	  DO xi_idx=1,NXIT,1  
            DFDZ(ns,nj_idx)=DFDZ(ns,nj_idx)+DFDXI(ns,xi_idx)*DXIDZ(xi_idx,nj_idx)
          ENDDO
        ENDDO
      ENDDO		  
    ENDDO	    
    
    DEALLOCATE(DFDXI)
    DEALLOCATE(DXIDZ)
    DEALLOCATE(DZDXI)    

    CALL EXITS("FINITE_ELASTICITY_GAUSS_DFDZ")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_DFDZ",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DFDZ")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,NEXT_NUMBER,NUMBER_OF_COMPONENTS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SUBTYPE)
        SELECT CASE(SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_SET%LINEARITY=EQUATIONS_SET_NONLINEAR
            EQUATIONS_SET%TIME_TYPE=EQUATIONS_SET_STATIC
            EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
!!TODO: maybe given negative user numbers to openCMISS generated fields???
            CALL FIELD_NEXT_NUMBER_FIND(EQUATIONS_SET%REGION,NEXT_NUMBER,ERR,ERROR,*999)
            CALL FIELD_CREATE_START(NEXT_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            CALL FIELD_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)
            CALL FIELD_DEPENDENT_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
            CALL FIELD_MESH_DECOMPOSITION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD% &
              & DECOMPOSITION,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
            NUMBER_OF_COMPONENTS=EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)% &
              & NUMBER_OF_COMPONENTS
            CALL FIELD_NUMBER_OF_COMPONENTS_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            !Default to the geometric interpolation setup
            DO component_idx=1,NUMBER_OF_COMPONENTS
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
                & component_idx,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)%COMPONENTS( &
                & component_idx)%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_NORMAL_VARIABLE_TYPE, &
                & component_idx,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)%COMPONENTS( &
                & component_idx)%MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
            ENDDO !component_idx
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO component_idx=1,NUMBER_OF_COMPONENTS
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
                  & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_NORMAL_VARIABLE_TYPE, &
                  & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
              ENDDO !component_idx
              CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD% &
                & SCALINGS%SCALING_TYPE,ERR,ERROR,*999)
            CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL FIELD_CREATE_FINISH(EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Do nothing
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FIXED_CONDITIONS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%FIXED_CONDITIONS)) THEN
              IF(EQUATIONS_SET%FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
                !Do nothing
                !?Initialise problem solution???
              ELSE
                CALL FLAG_ERROR("Equations set fixed conditions has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set fixed conditions is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,0,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_NORMAL_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_NO_SUBTYPE)        
        EQUATIONS_SET%CLASS=EQUATIONS_SET_ELASTICITY_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_FINITE_ELASTICITY_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SUBTYPE       
        CALL FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE, &
          & EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity problem.
  SUBROUTINE FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Laplace equation on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE)
        SELECT CASE(SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !DO nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            PROBLEM%NUMBER_OF_SOLUTIONS=1
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
!!TODO:
          CASE(PROBLEM_SETUP_FINISH_ACTION)
!!TODO:
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLUTION_TYPE)
          SOLUTION=>PROBLEM%SOLUTIONS(1)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SELECT CASE(ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              SOLUTION%LINEARITY=PROBLEM_SOLUTION_NONLINEAR
              CALL SOLUTION_MAPPING_CREATE_START(SOLUTION,SOLUTION_MAPPING,ERR,ERROR,*999)
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
              IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
                CALL SOLUTION_MAPPING_CREATE_FINISH(SOLUTION_MAPPING,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_SETUP_DO_ACTION)
              EQUATIONS_SET=>SOLUTION%EQUATIONS_SET_TO_ADD
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                !Check the equations set is from a finite elasticity equations set
                IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                  & EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SUBTYPE) THEN
                  SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING            
                  IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
                    CALL SOLUTION_MAPPING_EQUATIONS_SET_ADD(SOLUTION_MAPPING,SOLUTION%EQUATIONS_SET_TO_ADD,SOLUTION% &
                      & EQUATIONS_SET_ADDED_INDEX,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
               ELSE
                 CALL FLAG_ERROR("The equations set to add is not a finite elasticity equations set.",ERR,ERROR,*999)
               ENDIF
             ELSE
               CALL FLAG_ERROR("Equations set to add is not associated.",ERR,ERROR,*999)
             ENDIF
           CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity problem."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem solution is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_SETUP_SOLVER_TYPE)
          SOLUTION=>PROBLEM%SOLUTIONS(1)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SELECT CASE(ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              CALL SOLVER_CREATE_START(SOLUTION,SOLVER_NONLINEAR_TYPE,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              SOLVER=>SOLUTION%SOLVER
              IF(ASSOCIATED(SOLVER)) THEN                
                CALL SOLVER_CREATE_FINISH(SOLVER,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_SETUP_DO_ACTION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity problem."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem solution is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a finite elasticity type .
  SUBROUTINE FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_ELASTICITY_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE      
        CALL FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INITIAL_TYPE,PROBLEM_SETUP_START_ACTION, &
          & ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !
 
END MODULE FINITE_ELASTICITY_ROUTINES
