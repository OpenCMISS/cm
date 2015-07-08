!> \file
!> \author Chris Bradley
!> \brief This module contains all program wide constants.
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

!> This module contains all program wide constants.
MODULE CONSTANTS

  USE KINDS

  IMPLICIT NONE

  !Module parameters

  !> \addtogroup CONSTANTS_MathPhysicalConstants CONSTANTS::MathPhysicalConstants
  !> \see CONSTANTS
  !>@{ 
  REAL(DP), PARAMETER :: EULER=2.718281828459045235360287471352662497757_DP !<The double precision value of e \see CONSTANTS_MathPhysicalConstants,CONSTANTS
  REAL(DP), PARAMETER ::    PI=3.141592653589793238462643383279502884197_DP !<The double precision value of pi \see CONSTANTS_MathPhysicalConstants,CONSTANTS
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_DP !<The double value of 2pi \see CONSTANTS_MathPhysicalConstants,CONSTANTS
  !>@}
  
  !> \addtogroup CONSTANTS_NumericalConstants CONSTANTS::NumericalConstants
  !> \see CONSTANTS
  !>@{ 
  REAL(DP), PARAMETER :: CONVERGENCE_TOLERANCE=5.0_DP*EPSILON(1.0_DP) !<The convergence tolerance for double precision convergence calculations. Convergence tests should be of the form \f$\frac{|X_{i+1}-X_{i}|}{1+|X_{i}|}<\texttt{CONVERGENCE\_TOLERANCE}\f$ or for norms, \f$\frac{\|r\|}{\sqrt{n}+\|b\|}<\texttt{CONVERGENCE\_TOLERANCE}\f$
  !cpb 02/04/07 IBM compilers do not like this, initialise in BASE_ROUTINES_INITIALISE
  !REAL(DP), PARAMETER :: LOOSE_TOLERANCE=EPSILON(1.0_DP)**0.5
  REAL(DP) :: LOOSE_TOLERANCE !<The loose tolerance for double precision convergence calculations. Loose tolerance is to be used in the same manner as CONSTANTS::CONVERGENCE_TOLERANCE when a looser criterion is desired.
  REAL(DP), PARAMETER :: ZERO_TOLERANCE=5.0_DP*EPSILON(1.0_DP) !<The zero tolerance for double precision zero tests i.e., if(abs(x)>zero_tolerance) then...
  REAL(DP), PARAMETER :: CONVERGENCE_TOLERANCE_SP=5.0_SP*EPSILON(1.0_SP) !<The convergence tolerance for single precision convergence calculations. Convergence tests should be of the form \f$\frac{|X_{i+1}-X_{i}|}{1+|X_{i}|}<\texttt{CONVERGENCE\_TOLERANCE}\f$ or for norms, \f$\frac{\|r\|}{\sqrt{n}+\|b\|}<\texttt{CONVERGENCE\_TOLERANCE}\f$
  !cpb 02/04/07 IBM compilers do not like this, initialise in BASE_ROUTINES_INITIALISE
  !REAL(SP), PARAMETER :: LOOSE_TOLERANCE_SP=EPSILON(1.0_SP)**0.5
  REAL(SP) :: LOOSE_TOLERANCE_SP !<The loose tolerance for single precision convergence calculations. Loose tolerance is to be used in the same manner as CONSTANTS::CONVERGENCE_TOLERANCE_SP when a looser criterion is desired.
  REAL(SP), PARAMETER :: ZERO_TOLERANCE_SP=5.0_SP*EPSILON(1.0_SP) !<The zero tolerance for single precision zero tests i.e., if(abs(x)>zero_tolerance) then...
  !>@}

  !String parameters
  INTEGER(INTG), PARAMETER :: MAXSTRLEN=255 !<Maximum string length fro character strings

  !> \addtogroup CONSTANTS_DataTypes CONSTANTS::DataTypes
  !> Data type parameters for base data types
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: INTEGER_TYPE=1 !<Integer data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: SHORT_INTEGER_TYPE=2 !<Short integer data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LONG_INTEGER_TYPE=3 !<Long integer data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: SINGLE_REAL_TYPE=4 !<Single precision real data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: DOUBLE_REAL_TYPE=5 !<Double precision real data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: QUADRUPLE_REAL_TYPE=6 !<Quadruple precision real data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: CHARACTER_TYPE=7 !<Character data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LOGICAL_TYPE=8 !<Logical/boolean data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: SINGLE_COMPLEX_TYPE=9 !<Single precision complex data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: DOUBLE_COMPLEX_TYPE=10 !<Double precision complex data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: QUADRUPLE_COMPLEX_TYPE=11  !<Quadruple precision complex data type \see CONSTANTS_DataTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: C_INT_TYPE=12  !<C integer data type \see CONSTANTS_DataTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_EndianTypes CONSTANTS::EndianTypes
  !> Endian type parameters
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: BIG_ENDIAN_NUMBER=1 !<Big endian number type \see CONSTANTS_EndianTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LITTLE_ENDIAN_NUMBER=2 !<Little endian number type \see CONSTANTS_EndianTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_CharacterFormatTypes CONSTANTS::CharacterFormatTypes
  !> Bit format types for characters
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: ASCII_CHARACTER=1 !<ASCII character type \see CONSTANTS_CharacterFormatTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: UNICODE_CHARACTER=2 !<Unicode character type \see CONSTANTS_CharacterFormatTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_IntegerTypes CONSTANTS::IntegerTypes
  !> Bit format types for integers
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: TWOS_COMPLEMENT_INTEGER=1 !<Twos complement integer type \see CONSTANTS_IntegerFormatTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: SIGNED_MAGNITUDE_INTEGER=2 !<Signed magnitude integer type \see CONSTANTS_IntegerFormatTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_RealFormatTypes CONSTANTS::RealFormatTypes
  !> Bit format types for reals
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: SPIEEE_NUMBER=1 !<Single precision IEEE real type \see CONSTANTS_RealFormatTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: DPIEEE_NUMBER=2 !<Double precision IEEE real type \see CONSTANTS_RealFormatTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_ComputerSystemTypes CONSTANTS::ComputerSystemTypes
  !> Computer system type parameters
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: DEC_COMPUTER=1 !<Digital computer system type \see CONSTANTS_ComputerSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: SGI_COMPUTER=2 !<Silicon Graphics computer system type \see CONSTANTS_ComputerSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: IBM_COMPUTER=3 !<IBM system type \see CONSTANTS_ComputerSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: CRAY_COMPUTER=4 !<Cray computer system type \see CONSTANTS_ComputerSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: PC_COMPUTER=5 !<PC computer system type \see CONSTANTS_ComputerSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: UNKNOWN_COMPUTER=255 !<Unknown computer system type \see CONSTANTS_ComputerSystemTypes,CONSTANTS
  !>@}
  
  !> \addtogroup CONSTANTS_OperatingSystemTypes CONSTANTS::OperatingSystemTypes
  !> Operating system type parameters
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: VMS_OS=1 !<VMS operating system type \see CONSTANTS_OperatingSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: IRIX_OS=2 !<IRIX operating system type \see CONSTANTS_OperatingSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: WINDOWS_OS=3 !<Windows operating system type \see CONSTANTS_OperatingSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LINUX_OS=4 !<Linux operating system type \see CONSTANTS_OperatingSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: AIX_OS=5 !<AIX operating system type \see CONSTANTS_OperatingSystemTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: UNKNOWN_OS=255 !<Unknown operating system type \see CONSTANTS_OperatingSystemTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_LiblraryTypes CONSTANTS::LibraryTypes
  !> \brief Library type identifiers
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: LIBRARY_CMISS_TYPE=1 !<CMISS (internal) library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_PETSC_TYPE=2 !<PETSc library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_MUMPS_TYPE=3 !<MUMPS library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_SUPERLU_TYPE=4 !<SuperLU library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_SPOOLES_TYPE=5 !<SPOOLES library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_UMFPACK_TYPE=6 !<UMFPack library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_LUSOL_TYPE=7 !<LUSOL library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_ESSL_TYPE=8 !<ESSL library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_LAPACK_TYPE=9 !<LAPACK library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_TAO_TYPE=10 !<TAO library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_HYPRE_TYPE=11 !<Hypre library type \see CONSTANTS_LibraryTypes,CONSTANTS
  INTEGER(INTG), PARAMETER :: LIBRARY_PASTIX_TYPE=12 !<PaStiX library type \see CONSTANTS_LibraryTypes,CONSTANTS
  !>@}

  !> \addtogroup CONSTANTS_PartialDerivativeConstants CONSTANTS::PartialDerivativeConstants
  !> \brief Partial derivative constant identifiers
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: NO_PART_DERIV=1 !<No partial derivative i.e., u \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: FIRST_PART_DERIV=2 !<First partial derivative i.e., du/ds \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: SECOND_PART_DERIV=3 !<Second partial derivative i.e., d^2u/ds^2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: THIRD_PART_DERIV=4 !<Third partial derivative i.e., d^3u/ds^3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1=2 !<First partial derivative in the s1 direction i.e., du/ds1 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S1=3 !<Second partial derivative in the s1 direction i.e., d^2u/ds1ds1 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S2=4 !<First partial derivative in the s2 direction i.e., du/ds2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S2_S2=5 !<Second partial derivative in the s2 direction i.e., d^2u/ds2ds2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S2=6 !<Cross derivative in the s1 and s2 direction i.e., d^2u/ds1ds2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S3=7 !<First partial derivative in the s3 direction i.e., du/ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S3_S3=8 !<Second partial derivative in the s3 direction i.e., d^2u/ds3ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S3=9 !<Cross derivative in the s1 and s3 direction i.e., d^2u/ds1ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S2_S3=10 !<Cross derivative in the s2 and s3 direction i.e., d^2u/ds2ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S2_S3=11 !<Cross derivative in the s1, s2 and s3 direction i.e., d^3u/ds1ds2ds3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S4=12 !<First partial derivative in the s4 direction i.e., du/ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S4_S4=13 !<Second partial derivative in the s4 direction i.e., d^2u/ds4ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S4=14 !<Cross derivative in the s1 and s4 direction i.e., d^2u/ds1ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S2_S4=15 !<Cross derivative in the s2 and s4 direction i.e., d^2u/ds2ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S3_S4=16 !<Cross derivative in the s3 and s4 direction i.e., d^2u/ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S2_S4=17 !<Cross derivative in the s1, s2 and s4 direction i.e., d^3u/ds1ds2ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S3_S4=18 !<Cross derivative in the s1, s3 and s4 direction i.e., d^3u/ds1ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S2_S3_S4=19 !<Cross derivative in the s2, s3 and s4 direction i.e., d^3u/ds2ds3ds4 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: PART_DERIV_S1_S4_S4=20 !<Cross derivative in the s2, s4 and s4 direction i.e., d^3u/ds1ds4^2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  
  INTEGER(INTG), PARAMETER :: PART_DERIV_S2_S4_S4=21 !<Cross derivative in the s2, s4 and s4 direction i.e., d^3u/ds2ds4^2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  
  INTEGER(INTG), PARAMETER :: PART_DERIV_S3_S4_S4=22 !<Cross derivative in the s3, s4 and s4 direction i.e., d^3u/ds3ds4^2 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  
  INTEGER(INTG), PARAMETER :: PART_DERIV_S4_S4_S4=23 !<Third partial derivative in the s4 direction i.e., d^3u/ds4^3 \see CONSTANTS_PartialDerivativeConstants,CONSTANTS
  !>@}
  
  !> \addtogroup CONSTANTS_GlobalDerivativeConstants CONSTANTS::GlobalDerivativeConstants
  !> \brief Global derivative constant identifiers
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: MAXIMUM_GLOBAL_DERIV_NUMBER=8 !<The maximum global derivative number
  INTEGER(INTG), PARAMETER :: NO_GLOBAL_DERIV=1 !<No global derivative i.e., u \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S1=2 !<First global derivative in the s1 direction i.e., du/ds1 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S2=3 !<First global derivative in the s2 direction i.e., du/ds2 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S1_S2=4 !<Global Cross derivative in the s1 and s2 direction i.e., d^2u/ds1ds2 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S3=5 !<First global derivative in the s3 direction i.e., du/ds3 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S1_S3=6 !<Global Cross derivative in the s1 and s3 direction i.e., d^2u/ds1ds3 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S2_S3=7 !<Global Cross derivative in the s2 and s3 direction i.e., d^2u/ds2ds3 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GLOBAL_DERIV_S1_S2_S3=8 !<Cross derivative in the s1, s2 and s3 direction i.e., d^3u/ds1ds2ds3 \see CONSTANTS_GlobalDerivativeConstants,CONSTANTS
  !>@}
  
  !> \addtogroup CONSTANTS_PhysicalDerivativeConstants CONSTANTS::PhysicalDerivativeConstants
  !> \brief Physical derivative constant identifiers
  !> \see CONSTANTS
  !>@{ 
  INTEGER(INTG), PARAMETER :: MAXIMUM_PHYSICAL_DERIV_NUMBER=2 !<The maximum physical derivative number
  INTEGER(INTG), PARAMETER :: NO_PHYSICAL_DERIV=1 !<No physical derivative i.e., u \see CONSTANTS_PhysicalDerivativeConstants,CONSTANTS
  INTEGER(INTG), PARAMETER :: GRADIENT_PHYSICAL_DERIV=2 !<Gradient physical derivative i.e., grad u \see CONSTANTS_PhysicalDerivativeConstants,CONSTANTS
  !>@}
  
  INTEGER(INTG) :: PARTIAL_DERIVATIVE_INDEX(23,4) = RESHAPE( &
    & [ NO_PART_DERIV,FIRST_PART_DERIV,SECOND_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV, &
    &    NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV,SECOND_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV, &
    &    NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,FIRST_PART_DERIV,SECOND_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    FIRST_PART_DERIV,NO_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV,NO_PART_DERIV, &
    &    NO_PART_DERIV,FIRST_PART_DERIV,SECOND_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV, &
    &    FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,SECOND_PART_DERIV, &
    &    SECOND_PART_DERIV,SECOND_PART_DERIV,THIRD_PART_DERIV ], [23,4]) !<Partial derivative index map. PARTIAL_DERIVATIVE_INDEX(idx,nic) gives the order of the partial derivative in the ni(c)'th direction for the idx'th partial derivative value.

  INTEGER(INTG) :: PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(4) = [ PART_DERIV_S1,PART_DERIV_S2,PART_DERIV_S3,PART_DERIV_S4 ] !<PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(nic) gives the partial derivative index for the first derivative in the ni(c)'th direction
  
  INTEGER(INTG) :: PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(4) = [ PART_DERIV_S1_S1,PART_DERIV_S2_S2,PART_DERIV_S3_S3, &
    & PART_DERIV_S4_S4 ] !<PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(nic) gives the partial derivative index for the second derivative in the ni(c)'th direction

  INTEGER(INTG) :: PARTIAL_DERIVATIVE_MAXIMUM_MAP(4) = [ PART_DERIV_S1_S1,PART_DERIV_S1_S2,PART_DERIV_S1_S2_S3, &
    & PART_DERIV_S4_S4_S4 ] !<PARTIAL_DERIVATIVE_MAXIMUM_MAP(nic) gives the maximum of partial derivative index for the the ni(c)'th direction

  INTEGER(INTG) :: PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(20) = [ NO_GLOBAL_DERIV,GLOBAL_DERIV_S1,0,GLOBAL_DERIV_S2,0, &
    & GLOBAL_DERIV_S1_S2,GLOBAL_DERIV_S3,0,GLOBAL_DERIV_S1_S3,GLOBAL_DERIV_S2_S3,GLOBAL_DERIV_S1_S2_S3,0,0,0,0,0,0,0,0,0 ] !<PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(nu) gives the global derivative index for the the nu'th partial derivative. If no global derivative exists the map is zero

  INTEGER(INTG) :: GLOBAL_DERIVATIVE_PARTIAL_DERIVATIVE_MAP(8) = [ NO_PART_DERIV,PART_DERIV_S1,PART_DERIV_S2,PART_DERIV_S1_S2, &
    & PART_DERIV_S3,PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3] !<GLOBAL_DERIVATIVE_PARTIAL_DERIVATIVE_MAP(nk) gives the partial derivative index for the the nk'th global derivative.

  !>
  INTEGER(INTG) :: GLOBAL_DERIVATIVE_MAXIMUM_MAP(3) = [ GLOBAL_DERIV_S1,GLOBAL_DERIV_S1_S2,GLOBAL_DERIV_S1_S2_S3 ] !<GLOBAL_DERIVATIVE_MAXIMUM_MAP(ni) gives the maximum of global derivative index for the the ni'th direction
  
  !Other xi directions
  !>
  INTEGER(INTG) :: OTHER_XI_DIRECTIONS2(2) = [ 2,1 ] !<OTHER_XI_DIRECTIONS2(ni) gives the other xi direction for direction ni for a two dimensional element
  !>
  INTEGER(INTG) :: OTHER_XI_DIRECTIONS3(3,3,2) = RESHAPE([ 1,2,3,2,1,1,3,3,2,0,3,2,3,0,1,2,1,0 ],[3,3,2]) !<OTHER_XI_DIRECTIONS3(ni,nii,type) gives the other xi directions for direction ni for a three dimensional element. When type=1 then the nii index gives the other two xi directions (for nii=2,3) and when type=2 then ni and nii are used to give the third xi direction.

  INTEGER(INTG) :: OTHER_XI_DIRECTIONS4(4,3) = RESHAPE([ 2,3,4,1,3,4,1,2,4,1,2,3 ],[4,3]) !<OTHER_XI_DIRECTIONS4(nic,nii) gives the other xi coordinates for coordinate nic for a simplex element.
        
  INTEGER(INTG) :: OTHER_XI_ORIENTATIONS2(2) = [1,-1] !<OTHER_XI_ORIENTATIONSS2(ni) gives the orientation of the given xi direction and the other xi direction. Is equal to leviCivita(ni,OTHER_XI_DIRECTIONS2(ni)) where leviCivita is the Levi-Civita or alternating symbol
  INTEGER(INTG) :: OTHER_XI_ORIENTATIONS3(3,3) = RESHAPE([0,-1,1,1,0,-1,-1,1,0],[3,3]) !<OTHER_XI_ORIENTATIONSS3(ni,nii) gives the orientation of the given two xi directions. Is equal to leviCivita(ni,nii,OTHER_XI_DIRECTIONS3(ni,nii,2)) where leviCivita is the Levi-Civita or alternating symbol

  !Tensor to voigt indices and voigt to tensor indices maps.
  !>
  INTEGER(INTG), PARAMETER :: TENSOR_TO_VOIGT2(2,2)=RESHAPE([1,3,3,2],[2,2]) !TENSOR_TO_VOIGT2(i,j) converts a pair (i,j) of rank 2 symmetric tensor indices to Voigt index (a) in 2 dimensions.
  INTEGER(INTG), PARAMETER :: VOIGT_TO_TENSOR2(2,3)=RESHAPE([1,1,2,2,1,2],[2,3]) !VOIGT_TO_TENSOR2(a) converts a Voigt index (a) to a pair (i,j) of rank 2 symmetric tensor indices in 2 dimensions.
  INTEGER(INTG), PARAMETER :: TENSOR_TO_VOIGT3(3,3)=RESHAPE([1,4,5,4,2,6,5,6,3],[3,3]) !TENSOR_TO_VOIGT3(i,j) converts a pair (i,j) of rank 2 symmetric tensor indices to Voigt index (a) in 3 dimensions.
  INTEGER(INTG), PARAMETER :: VOIGT_TO_TENSOR3(2,6)=RESHAPE([1,1,2,2,3,3,1,2,1,3,2,3],[2,6]) !VOIGT_TO_TENSOR3(i,j) converts a Voigt index (a) to a pair (i,j) of rank 2 symmetric tensor indices in 3 dimensions.

END MODULE CONSTANTS
