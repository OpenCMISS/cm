!> \file
!> \author Chris Bradley
!> \brief This module contains all routines dealing with (non-distributed) matrix and vectors types.
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
!> \section MATRIX_VECTOR_MatrixStorageStructures MATRIX STORAGE STRUCTURES
!>    The matrix storage structures used are governed by the STORAGE parameter associated with the array. If STORAGE is
!>    MATRIX_BLOCK_STORAGE_TYPE the matrix is not sparse and the the non-sparse matrix dimension M is used to calculate
!>    the matrix storage locations. If storage is MATRIX_DIAGONAL_STORAGE_TYPE then only the matrix diagonal is stored.
!>    If storage is MATRIX_COLUMN_MAJOR_STORAGE_TYPE the matrix is not sparse and the non-sparse
!>    matrix dimension MAX_M (>=M) is used to calcualte the matrix storage locations. If storage is
!>    MATRIX_ROW_MAJOR_STORAGE_TYPE the matrix is not sparse and the non-sparse matrix dimension MAX_N (>=N) is used to
!>    calcualte the matrix storage locations. If STORAGE is MATRIX_COMPRESSED_ROW_STORAGE_TYPE the matrix has compressed row
!>    storage/sparsity (see below) and the sparsity structure arrays ROW_INDICES and COLUMN_INDICES are used for the
!>    storage location calculation. If STORAGE is MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE the matrix has compressed column
!>    storage/sparsity (see below) and the sparsity structure arrays ROW_INDICES and COLUMN_INDICES are used for the
!>    storage location calculation. If STORAGE is MATRIX_ROW_COLUMN_STORAGE_TYPE the matrix has row column
!>    storage/sparsity (see below) and the sparsity structure arrays ROW_INDICES and COLUMN_INDICES are used for the
!>    storage location calculation. 
!>    
!>    \subsection MATRIX_VECTOR_CompressedRowStorage COMPRESSED-ROW STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a one dimensional array of length SIZE
!>      (=NUMBER_NON_ZEROS) (where NUMBER_NON_ZEROS=sxMxN, s is the sparsity of the array) that stores only the non-zero
!>      elements of the matrix. Two additional arrays ROW_INDICES and COLUMN_INDICES store the positions of the non-zero
!>      elements. ROW_INDICES is of length M+1 and COLUMN is of length NUMBER_NON_ZEROS. ROW_INDICES(i) stores the position
!>      in COLUMN_INDICES of  the start of row i. The M+1 position of ROW_INDICES stores the size of COLUMN_INDICES+1 i.e.,
!>      NUMBER_NON_ZEROS+1. The number of non-zero elements in row i can be found from ROW_INDICES(i+1)-ROW_INDICES(i).
!>      COLUMN_INDICES(nz) gives the column number for non-zero element nz. See also COMPRESSED-COLUMN storage.
!>   
!>      Example of the compressed-row storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%
!>      \verbatim
!>      
!>      GX  1 2 3 4 5 6        
!>         ____________        DATA(nz)                        
!>       1| 0 A 0 B 0 0          A B C D E F G H  
!>       2| 0 0 C 0 0 0          
!>       3| 0 0 0 0 D E        ROW_INDICES(i)
!>       4| F 0 0 0 0 0          1 3 4 6 7 8 9
!>       5| 0 0 G 0 0 0        COLUMN_INDICES(i)  
!>       6| 0 0 0 0 0 H          2 4 3 5 6 1 3 6
!>      
!>      \endverbatim
!>
!>    \subsection MATRIX_VECTOR_CompressedColumnStorage COMPRESSED-COLUMN STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a one dimensional array of length SIZE
!>      (=NUMBER_NON_ZEROS) (where NUMBER_NON_ZEROS=sxMxN, s is the sparsity of the array) that stores only the non-zero
!>      elements of the matrix. Two additional arrays ROW_INDICES and COLUMN_INDICES store the positions of the non-zero
!>      elements. ROW_INDICES is of length NUMBER_NON_ZEROS and COLUMN is of length N+1. COLUMN_INDICES(j) stores the position
!>      in ROW_INDICES of  the start of column j. The N+1 position of COLUMN_INDICES stores the size of ROW_INDICES+1 i.e.,
!>      NUMBER_NON_ZEROS+1. The number of non-zero elements in column j can be found from COLUMN_INDICES(j+1)-COLUMN_INDICES(j).
!>      ROW_INDICES(nz) gives the row number for non-zero element nz. See also COMPRESSED-ROW storage.
!>   
!>      Example of compressed-column storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%
!>      \verbatim
!>      
!>      GX  1 2 3 4 5 6        
!>         ____________        DATA(nz)                        
!>       1| 0 A 0 B 0 0          F A C G B D E H  
!>       2| 0 0 C 0 0 0          
!>       3| 0 0 0 0 D E        ROW_INDICES(i)
!>       4| F 0 0 0 0 0          4 1 2 5 1 3 3 6
!>       5| 0 0 G 0 0 0        COLUMN_INDICES(i)  
!>       6| 0 0 0 0 0 H          1 2 3 5 6 7 9
!>      
!>      \endverbatim
!>
!>    \subsection MATRIX_VECTOR_RowColumnStorage ROW-COLUMN STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a one dimensional array of length SIZE
!>      (=NUMBER_NON_ZEROS) (where NUMBER_NON_ZEROS=sxMxN, s is the sparsity of the array) that stores only the non-zero
!>      elements of the matrix. Two additional arrays ROW_INDICES and COLUMN_INDICES store the positions of the non-zero
!>      elements. Both ROW_INDICES and COLUMN_INDICES are of length NUMBER_NON_ZEROS. ROW_INDICES(nz) gives the row number for
!>      non-zero element nz and COLUMN_INDICES(nz) gives the column number for non-zero element nz.
!>  
!>      Example of row-column storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%
!>      \verbatim
!>   
!>      GX  1 2 3 4 5 6        
!>         ____________        DATA(nz)                        
!>       1| 0 A 0 B 0 0          A B C D E F G H  
!>       2| 0 0 C 0 0 0          
!>       3| 0 0 0 0 D E        ROW_INDICES(i)
!>       4| F 0 0 0 0 0          1 1 2 3 3 4 5 6 
!>       5| 0 0 G 0 0 0        COLUMN_INDICES(i)  
!>       6| 0 0 0 0 0 H          2 4 3 5 6 1 3 6
!>
!>      \endverbatim

!>This module contains all routines dealing with (non-distributed) matrix and vectors types.
MODULE MATRIX_VECTOR

  USE BASE_ROUTINES
  USE CONSTANTS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE STRINGS
  USE TYPES
  USE LINKEDLIST_ROUTINES

  IMPLICIT NONE

  !PRIVATE

  !Module parameters

  !> \addtogroup MATRIX_VECTOR_DataTypes MATRIX_VECTOR::DataTypes
  !> \brief Matrix vector data types
  !> \see MATRIX_VECTOR
  !>@{
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_INTG_TYPE=INTEGER_TYPE !<Integer matrix-vector data type \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_SP_TYPE=SINGLE_REAL_TYPE !<Single precision real matrix-vector data type \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_DP_TYPE=DOUBLE_REAL_TYPE !<Double precision real matrix-vector data type \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_L_TYPE=LOGICAL_TYPE !<Logical matrix-vector data type \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
  !>@}
  
  !> \addtogroup MATRIX_VECTOR_StorageTypes MATRIX_VECTOR::StorageTypes
  !> \brief Matrix-vector storage type parameters
  !> \see MATRIX_VECTOR_MatrixStorageStructures,MATRIX_VECTOR
  !>@{
  INTEGER(INTG), PARAMETER :: MATRIX_BLOCK_STORAGE_TYPE=0 !<Matrix block storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_DIAGONAL_STORAGE_TYPE=1 !<Matrix diagonal storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_COLUMN_MAJOR_STORAGE_TYPE=2 !<Matrix column major storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_ROW_MAJOR_STORAGE_TYPE=3 !<Matrix row major storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_COMPRESSED_ROW_STORAGE_TYPE=4 !<Matrix compressed row storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE=5 !<Matrix compressed column storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: MATRIX_ROW_COLUMN_STORAGE_TYPE=6 !<Matrix row-column storage type \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
  !>@}
  
  !Module types

  !Matrix types
  
  !Module variables

  INTEGER(INTG), SAVE :: MATRIX_VECTOR_ID=1

  !Interfaces

  INTERFACE MATRIX_ALL_VALUES_SET
    MODULE PROCEDURE MATRIX_ALL_VALUES_SET_INTG
    MODULE PROCEDURE MATRIX_ALL_VALUES_SET_SP
    MODULE PROCEDURE MATRIX_ALL_VALUES_SET_DP
    MODULE PROCEDURE MATRIX_ALL_VALUES_SET_L
  END INTERFACE !MATRIX_ALL_VALUES_SET

  INTERFACE MATRIX_DATA_GET
    MODULE PROCEDURE MATRIX_DATA_GET_INTG
    MODULE PROCEDURE MATRIX_DATA_GET_SP
    MODULE PROCEDURE MATRIX_DATA_GET_DP
    MODULE PROCEDURE MATRIX_DATA_GET_L
  END INTERFACE !MATRIX_DATA_GET

  INTERFACE MATRIX_VALUES_ADD
    MODULE PROCEDURE MATRIX_VALUES_ADD_INTG
    MODULE PROCEDURE MATRIX_VALUES_ADD_INTG1
    MODULE PROCEDURE MATRIX_VALUES_ADD_INTG2
    MODULE PROCEDURE MATRIX_VALUES_ADD_SP
    MODULE PROCEDURE MATRIX_VALUES_ADD_SP1
    MODULE PROCEDURE MATRIX_VALUES_ADD_SP2
    MODULE PROCEDURE MATRIX_VALUES_ADD_DP
    MODULE PROCEDURE MATRIX_VALUES_ADD_DP1
    MODULE PROCEDURE MATRIX_VALUES_ADD_DP2
    MODULE PROCEDURE MATRIX_VALUES_ADD_L
    MODULE PROCEDURE MATRIX_VALUES_ADD_L1
    MODULE PROCEDURE MATRIX_VALUES_ADD_L2
  END INTERFACE !MATRIX_VALUES_ADD

  INTERFACE MATRIX_VALUES_GET
    MODULE PROCEDURE MATRIX_VALUES_GET_INTG
    MODULE PROCEDURE MATRIX_VALUES_GET_INTG1
    MODULE PROCEDURE MATRIX_VALUES_GET_INTG2
    MODULE PROCEDURE MATRIX_VALUES_GET_SP
    MODULE PROCEDURE MATRIX_VALUES_GET_SP1
    MODULE PROCEDURE MATRIX_VALUES_GET_SP2
    MODULE PROCEDURE MATRIX_VALUES_GET_DP
    MODULE PROCEDURE MATRIX_VALUES_GET_DP1
    MODULE PROCEDURE MATRIX_VALUES_GET_DP2
    MODULE PROCEDURE MATRIX_VALUES_GET_L
    MODULE PROCEDURE MATRIX_VALUES_GET_L1
    MODULE PROCEDURE MATRIX_VALUES_GET_L2
  END INTERFACE !MATRIX_VALUES_GET
  
  INTERFACE MATRIX_VALUES_SET
    MODULE PROCEDURE MATRIX_VALUES_SET_INTG
    MODULE PROCEDURE MATRIX_VALUES_SET_INTG1
    MODULE PROCEDURE MATRIX_VALUES_SET_INTG2
    MODULE PROCEDURE MATRIX_VALUES_SET_SP
    MODULE PROCEDURE MATRIX_VALUES_SET_SP1
    MODULE PROCEDURE MATRIX_VALUES_SET_SP2
    MODULE PROCEDURE MATRIX_VALUES_SET_DP
    MODULE PROCEDURE MATRIX_VALUES_SET_DP1
    MODULE PROCEDURE MATRIX_VALUES_SET_DP2
    MODULE PROCEDURE MATRIX_VALUES_SET_L
    MODULE PROCEDURE MATRIX_VALUES_SET_L1
    MODULE PROCEDURE MATRIX_VALUES_SET_L2
  END INTERFACE !MATRIX_VALUES_SET

  INTERFACE VECTOR_ALL_VALUES_SET
    MODULE PROCEDURE VECTOR_ALL_VALUES_SET_INTG
    MODULE PROCEDURE VECTOR_ALL_VALUES_SET_SP
    MODULE PROCEDURE VECTOR_ALL_VALUES_SET_DP
    MODULE PROCEDURE VECTOR_ALL_VALUES_SET_L
  END INTERFACE !VECTOR_ALL_VALUES_SET

  INTERFACE VECTOR_DATA_GET
    MODULE PROCEDURE VECTOR_DATA_GET_INTG
    MODULE PROCEDURE VECTOR_DATA_GET_SP
    MODULE PROCEDURE VECTOR_DATA_GET_DP
    MODULE PROCEDURE VECTOR_DATA_GET_L
  END INTERFACE !VECTOR_DATA_GET

   INTERFACE VECTOR_VALUES_GET
    MODULE PROCEDURE VECTOR_VALUES_GET_INTG
    MODULE PROCEDURE VECTOR_VALUES_GET_INTG1
    MODULE PROCEDURE VECTOR_VALUES_GET_SP
    MODULE PROCEDURE VECTOR_VALUES_GET_SP1
    MODULE PROCEDURE VECTOR_VALUES_GET_DP
    MODULE PROCEDURE VECTOR_VALUES_GET_DP1
    MODULE PROCEDURE VECTOR_VALUES_GET_L
    MODULE PROCEDURE VECTOR_VALUES_GET_L1
  END INTERFACE !VECTOR_VALUES_GET
  
  INTERFACE VECTOR_VALUES_SET
    MODULE PROCEDURE VECTOR_VALUES_SET_INTG
    MODULE PROCEDURE VECTOR_VALUES_SET_INTG1
    MODULE PROCEDURE VECTOR_VALUES_SET_SP
    MODULE PROCEDURE VECTOR_VALUES_SET_SP1
    MODULE PROCEDURE VECTOR_VALUES_SET_DP
    MODULE PROCEDURE VECTOR_VALUES_SET_DP1
    MODULE PROCEDURE VECTOR_VALUES_SET_L
    MODULE PROCEDURE VECTOR_VALUES_SET_L1
  END INTERFACE !VECTOR_VALUES_SET

  PUBLIC MATRIX_VECTOR_INTG_TYPE,MATRIX_VECTOR_SP_TYPE,MATRIX_VECTOR_DP_TYPE,MATRIX_VECTOR_L_TYPE

  PUBLIC MATRIX_BLOCK_STORAGE_TYPE,MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE, &
    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE

  PUBLIC MATRIX_CREATE_FINISH,MATRIX_CREATE_START,MATRIX_DATA_GET,MATRIX_DATA_TYPE_SET,MATRIX_DESTROY, &
    & MATRIX_DUPLICATE,MATRIX_MAX_COLUMNS_PER_ROW_GET,MATRIX_NUMBER_NON_ZEROS_SET,MATRIX_NUMBER_NON_ZEROS_GET,MATRIX_MAX_SIZE_SET, &
    & MATRIX_OUTPUT,MATRIX_SIZE_SET,MATRIX_STORAGE_LOCATION_FIND,MATRIX_STORAGE_LOCATIONS_SET,MATRIX_STORAGE_TYPE_GET, &
    & MATRIX_STORAGE_TYPE_SET,MATRIX_VALUES_ADD,MATRIX_VALUES_GET,MATRIX_VALUES_SET

  PUBLIC VECTOR_ALL_VALUES_SET,VECTOR_CREATE_FINISH,VECTOR_CREATE_START,VECTOR_DATA_GET,VECTOR_DATA_TYPE_SET,VECTOR_DESTROY, &
    & VECTOR_DUPLICATE,VECTOR_SIZE_SET,VECTOR_VALUES_GET,VECTOR_VALUES_SET

  PUBLIC MATRIX_LINKLIST_SET,MATRIX_LINKLIST_GET
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Sets all values in an integer matrix to the specified value.
  SUBROUTINE MATRIX_ALL_VALUES_SET_INTG(MATRIX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_ALL_VALUES_SET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          MATRIX%DATA_INTG=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_ALL_VALUES_SET_INTG")
    RETURN
999 CALL ERRORS("MATRIX_ALL_VALUES_SET_INTG",ERR,ERROR)
    CALL EXITS("MATRIX_ALL_VALUES_SET_INTG")
    RETURN 1
  END SUBROUTINE MATRIX_ALL_VALUES_SET_INTG

  !
  !================================================================================================================================
  !

  !>Sets all values in a single precision matrix to the specified value.
  SUBROUTINE MATRIX_ALL_VALUES_SET_SP(MATRIX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    REAL(SP), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_ALL_VALUES_SET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          MATRIX%DATA_SP=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_ALL_VALUES_SET_SP")
    RETURN
999 CALL ERRORS("MATRIX_ALL_VALUES_SET_SP",ERR,ERROR)
    CALL EXITS("MATRIX_ALL_VALUES_SET_SP")
    RETURN 1
  END SUBROUTINE MATRIX_ALL_VALUES_SET_SP

  !
  !================================================================================================================================
  !

  !>Sets all values in a double precision matrix to the specified value.
  SUBROUTINE MATRIX_ALL_VALUES_SET_DP(MATRIX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_ALL_VALUES_SET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          MATRIX%DATA_DP=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_ALL_VALUES_SET_DP")
    RETURN
999 CALL ERRORS("MATRIX_ALL_VALUES_SET_DP",ERR,ERROR)
    CALL EXITS("MATRIX_ALL_VALUES_SET_DP")
    RETURN 1
  END SUBROUTINE MATRIX_ALL_VALUES_SET_DP

  !
  !================================================================================================================================
  !

  !>Sets all values in a logical matrix to the specified value.
  SUBROUTINE MATRIX_ALL_VALUES_SET_L(MATRIX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    LOGICAL, INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_ALL_VALUES_SET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          MATRIX%DATA_L=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_ALL_VALUES_SET_L")
    RETURN
999 CALL ERRORS("MATRIX_ALL_VALUES_SET_L",ERR,ERROR)
    CALL EXITS("MATRIX_ALL_VALUES_SET_L")
    RETURN 1
  END SUBROUTINE MATRIX_ALL_VALUES_SET_L

  !
  !================================================================================================================================
  !

  !>Finishes the creation a matrix.
  SUBROUTINE MATRIX_CREATE_FINISH(MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to finish
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,COUNT,row_idx,row_idx2
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("Matrix has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          IF(MATRIX%MAX_M==-1) MATRIX%MAX_M=MATRIX%M
          IF(MATRIX%MAX_N==-1) MATRIX%MAX_N=MATRIX%N
          MATRIX%SIZE=MATRIX%M*MATRIX%N
          MATRIX%NUMBER_NON_ZEROS=MATRIX%M*MATRIX%N
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=MATRIX%N
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          IF(MATRIX%MAX_M==-1) MATRIX%MAX_M=MATRIX%M
          IF(MATRIX%MAX_N==-1) MATRIX%MAX_N=MATRIX%N
          MATRIX%SIZE=MATRIX%M
          MATRIX%NUMBER_NON_ZEROS=MATRIX%M
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=1
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          IF(MATRIX%MAX_M==-1) CALL FLAG_ERROR("Maximum number of rows has not been set for this matrix.",ERR,ERROR,*999)
          IF(MATRIX%MAX_N==-1) CALL FLAG_ERROR("Maximum number of columns has not been set for this matrix.",ERR,ERROR,*999)
          MATRIX%SIZE=MATRIX%MAX_M*MATRIX%N
          MATRIX%NUMBER_NON_ZEROS=MATRIX%M*MATRIX%N
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=MATRIX%N
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          IF(MATRIX%MAX_M==-1) CALL FLAG_ERROR("Maximum number of rows has not been set for this matrix.",ERR,ERROR,*999)
          IF(MATRIX%MAX_N==-1) CALL FLAG_ERROR("Maximum number of columns has not been set for this matrix.",ERR,ERROR,*999)
          MATRIX%SIZE=MATRIX%M*MATRIX%MAX_N
          MATRIX%NUMBER_NON_ZEROS=MATRIX%M*MATRIX%N
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=MATRIX%N
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          IF(MATRIX%NUMBER_NON_ZEROS==-1) CALL FLAG_ERROR("Number of non-zeros has not been set for this matrix.",ERR,ERROR,*999)
          IF(MATRIX%MAX_M==-1) MATRIX%MAX_M=MATRIX%M
          IF(MATRIX%MAX_N==-1) MATRIX%MAX_N=MATRIX%N
          MATRIX%SIZE=MATRIX%NUMBER_NON_ZEROS
          IF(.NOT.ALLOCATED(MATRIX%COLUMN_INDICES))  &
            & CALL FLAG_ERROR("Matrix storage locations column indices have not been set.",ERR,ERROR,*999)
          IF(.NOT.ALLOCATED(MATRIX%ROW_INDICES))  &
            & CALL FLAG_ERROR("Matrix storage locations row indices have not been set.",ERR,ERROR,*999)
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=0
          DO row_idx=1,MATRIX%M
            IF((MATRIX%ROW_INDICES(row_idx+1)-MATRIX%ROW_INDICES(row_idx))>MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW) &
              & MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=MATRIX%ROW_INDICES(row_idx+1)-MATRIX%ROW_INDICES(row_idx)
          ENDDO !row_idx
        CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          IF(MATRIX%NUMBER_NON_ZEROS==-1) CALL FLAG_ERROR("Number of non-zeros has not been set for this matrix.",ERR,ERROR,*999)
          IF(MATRIX%MAX_M==-1) MATRIX%MAX_M=MATRIX%M
          IF(MATRIX%MAX_N==-1) MATRIX%MAX_N=MATRIX%N
          MATRIX%SIZE=MATRIX%NUMBER_NON_ZEROS
          IF(.NOT.ALLOCATED(MATRIX%COLUMN_INDICES))  &
            & CALL FLAG_ERROR("Matrix storage locations column indices have not been set.",ERR,ERROR,*999)
          IF(.NOT.ALLOCATED(MATRIX%ROW_INDICES))  &
            & CALL FLAG_ERROR("Matrix storage locations row indices have not been set.",ERR,ERROR,*999)
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=0
          DO row_idx=1,MATRIX%M
            COUNT=0
            DO column_idx=1,MATRIX%N
              DO row_idx2=MATRIX%COLUMN_INDICES(column_idx),MATRIX%COLUMN_INDICES(column_idx+1)-1
                IF(MATRIX%ROW_INDICES(row_idx2)==row_idx) COUNT=COUNT+1
              ENDDO !row_idx2
            ENDDO !column_idx
            IF(COUNT>MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW) MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=COUNT
          ENDDO !row_idx
        CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
          IF(MATRIX%NUMBER_NON_ZEROS==-1) CALL FLAG_ERROR("Number of non-zeros has not been set for this matrix.",ERR,ERROR,*999)
          IF(MATRIX%MAX_M==-1) MATRIX%MAX_M=MATRIX%M
          IF(MATRIX%MAX_N==-1) MATRIX%MAX_N=MATRIX%N
          MATRIX%SIZE=MATRIX%NUMBER_NON_ZEROS  
          IF(.NOT.ALLOCATED(MATRIX%COLUMN_INDICES))  &
            & CALL FLAG_ERROR("Matrix storage locations column indices have not been set.",ERR,ERROR,*999)
          IF(.NOT.ALLOCATED(MATRIX%ROW_INDICES))  &
            & CALL FLAG_ERROR("Matrix storage locations row indices have not been set.",ERR,ERROR,*999)
          MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=0
          DO row_idx=1,MATRIX%M
            COUNT=0
            DO row_idx2=1,MATRIX%NUMBER_NON_ZEROS
              IF(MATRIX%ROW_INDICES(row_idx2)==row_idx) COUNT=COUNT+1
            ENDDO !row_idx2            
            IF(COUNT>MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW) MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=COUNT
          ENDDO !row_idx          
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(MATRIX%SIZE>0) THEN
          SELECT CASE(MATRIX%DATA_TYPE)
          CASE(MATRIX_VECTOR_INTG_TYPE)
            ALLOCATE(MATRIX%DATA_INTG(MATRIX%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix integer data.",ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_SP_TYPE)
            ALLOCATE(MATRIX%DATA_SP(MATRIX%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix single precision data.",ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_DP_TYPE)
            ALLOCATE(MATRIX%DATA_DP(MATRIX%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix double precision data.",ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_L_TYPE)
            ALLOCATE(MATRIX%DATA_L(MATRIX%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix logical data.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The matrix data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
        MATRIX%ID=MATRIX_VECTOR_ID
        MATRIX_VECTOR_ID=MATRIX_VECTOR_ID+1
        MATRIX%MATRIX_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_CREATE_FINISH")
    RETURN
!!TODO: deallocate on error
999 CALL ERRORS("MATRIX_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("MATRIX_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE MATRIX_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation a matrix.
  SUBROUTINE MATRIX_CREATE_START(MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MATRIX_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      CALL FLAG_ERROR("Matrix is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(MATRIX,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate the matrix.",ERR,ERROR,*999)
      CALL MATRIX_INITIALISE(MATRIX,ERR,ERROR,*999)
      !Set the defaults
      MATRIX%DATA_TYPE=MATRIX_VECTOR_DP_TYPE
      MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
    ENDIF
    
    CALL EXITS("MATRIX_CREATE_START")
    RETURN
999 IF(ASSOCIATED(MATRIX)) CALL MATRIX_FINALISE(MATRIX,ERR,ERROR,*998)
998 CALL ERRORS("MATRIX_CREATE_START",ERR,ERROR)
    CALL EXITS("MATRIX_CREATE_START")
    RETURN 1
  END SUBROUTINE MATRIX_CREATE_START

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of an integer matrix. Note: the values can be used for read operations but a MATRIX_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE MATRIX_DATA_GET_INTG(MATRIX,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), POINTER :: DATA(:) !<On return a pointer to the matrix data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_DATA_GET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)    
        IF(MATRIX%MATRIX_FINISHED) THEN
          IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
            DATA=>MATRIX%DATA_INTG
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the integer data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_DATA_GET_INTG")
    RETURN
999 CALL ERRORS("MATRIX_DATA_GET_INTG",ERR,ERROR)
    CALL EXITS("MATRIX_DATA_GET_INTG")
    RETURN 1
  END SUBROUTINE MATRIX_DATA_GET_INTG

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a single precision matrix. Note: the values can be used for read operations but aMATRIX_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE MATRIX_DATA_GET_SP(MATRIX,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    REAL(SP), POINTER :: DATA(:) !<On return a pointer to the matrix data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_DATA_GET_SP",ERR,ERROR,*999)

     IF(ASSOCIATED(MATRIX)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(MATRIX%MATRIX_FINISHED) THEN
          IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
            DATA=>MATRIX%DATA_SP
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the single precision data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_DATA_GET_SP")
    RETURN
999 CALL ERRORS("MATRIX_DATA_GET_SP",ERR,ERROR)
    CALL EXITS("MATRIX_DATA_GET_SP")
    RETURN 1
  END SUBROUTINE MATRIX_DATA_GET_SP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a double precision matrix. Note: the values can be used for read operations but a MATRIX_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE MATRIX_DATA_GET_DP(MATRIX,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    REAL(DP), POINTER :: DATA(:) !<On return a pointer to the matrix data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_DATA_GET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(MATRIX%MATRIX_FINISHED) THEN
          IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
            DATA=>MATRIX%DATA_DP
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the double precision data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_DATA_GET_DP")
    RETURN
999 CALL ERRORS("MATRIX_DATA_GET_DP",ERR,ERROR)
    CALL EXITS("MATRIX_DATA_GET_DP")
    RETURN 1
  END SUBROUTINE MATRIX_DATA_GET_DP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a logical matrix. Note: the values can be used for read operations but a MATRIX_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE MATRIX_DATA_GET_L(MATRIX,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    LOGICAL, POINTER :: DATA(:) !<On return a pointer to the matrix data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_DATA_GET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(MATRIX%MATRIX_FINISHED) THEN
          IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
            DATA=>MATRIX%DATA_L
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the logical data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_DATA_GET_L")
    RETURN
999 CALL ERRORS("MATRIX_DATA_GET_L",ERR,ERROR)
    CALL EXITS("MATRIX_DATA_GET_L")
    RETURN 1
  END SUBROUTINE MATRIX_DATA_GET_L

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type of a matrix.
  SUBROUTINE MATRIX_DATA_TYPE_SET(MATRIX,DATA_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix 
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type to set for the matrix. \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_DATA_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN      
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("The matrix has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DATA_TYPE)
        CASE(MATRIX_VECTOR_INTG_TYPE)
          MATRIX%DATA_TYPE=MATRIX_VECTOR_INTG_TYPE
        CASE(MATRIX_VECTOR_SP_TYPE)
          MATRIX%DATA_TYPE=MATRIX_VECTOR_SP_TYPE
        CASE(MATRIX_VECTOR_DP_TYPE)
          MATRIX%DATA_TYPE=MATRIX_VECTOR_DP_TYPE
        CASE(MATRIX_VECTOR_L_TYPE)
          MATRIX%DATA_TYPE=MATRIX_VECTOR_L_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The matrix vector data type of "//TRIM(NUMBER_TO_VSTRING(DATA_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_DATA_TYPE_SET")
    RETURN
999 CALL ERRORS("MATRIX_DATA_TYPE_SET",ERR,ERROR)
    CALL EXITS("MATRIX_DATA_TYPE_SET")
    RETURN 1
  END SUBROUTINE MATRIX_DATA_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Destroys a matrix 
  SUBROUTINE MATRIX_DESTROY(MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("MATRIX_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      CALL MATRIX_FINALISE(MATRIX,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_DESTROY")
    RETURN
999 CALL ERRORS("MATRIX_DESTROY",ERR,ERROR)
    CALL EXITS("MATRIX_DESTROY")
    RETURN 1
  END SUBROUTINE MATRIX_DESTROY

  !
  !================================================================================================================================
  !

  !>Duplicates the matrix and returns a pointer to the duplicated matrix in NEWMATRIX.
  SUBROUTINE MATRIX_DUPLICATE(MATRIX,NEW_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to duplicate
    TYPE(MATRIX_TYPE), POINTER :: NEW_MATRIX !<On return a pointer to a new duplicated matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MATRIX_DUPLICATE",ERR,ERROR,*998)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(ASSOCIATED(NEW_MATRIX)) THEN
        CALL FLAG_ERROR("New matrix is already associated.",ERR,ERROR,*998)
      ELSE
        CALL MATRIX_CREATE_START(NEW_MATRIX,ERR,ERROR,*999)
        CALL MATRIX_DATA_TYPE_SET(NEW_MATRIX,MATRIX%DATA_TYPE,ERR,ERROR,*999)
        CALL MATRIX_SIZE_SET(NEW_MATRIX,MATRIX%M,MATRIX%N,ERR,ERROR,*999)
        CALL MATRIX_STORAGE_TYPE_SET(NEW_MATRIX,MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_DIAGONAL_STORAGE_TYPE)
          !Do nothing
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL MATRIX_MAX_SIZE_SET(NEW_MATRIX,MATRIX%MAX_M,MATRIX%MAX_N,ERR,ERROR,*999)          
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL MATRIX_NUMBER_NON_ZEROS_SET(NEW_MATRIX,MATRIX%NUMBER_NON_ZEROS,ERR,ERROR,*999)
          CALL MATRIX_STORAGE_LOCATIONS_SET(NEW_MATRIX,MATRIX%ROW_INDICES,MATRIX%COLUMN_INDICES,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        CALL MATRIX_CREATE_FINISH(NEW_MATRIX,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("MATRIX_DUPLICATE")
    RETURN
999 CALL MATRIX_FINALISE(NEW_MATRIX,ERR,ERROR,*998)
998 CALL ERRORS("MATRIX_DUPLICATE",ERR,ERROR)
    CALL EXITS("MATRIX_DUPLICATE")
    RETURN 1
  END SUBROUTINE MATRIX_DUPLICATE

  !
  !================================================================================================================================
  !

  !>Finalises a matrix and deallocates all memory.
  SUBROUTINE MATRIX_FINALISE(MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(ALLOCATED(MATRIX%ROW_INDICES)) DEALLOCATE(MATRIX%ROW_INDICES)
      IF(ALLOCATED(MATRIX%COLUMN_INDICES)) DEALLOCATE(MATRIX%COLUMN_INDICES)
      IF(ALLOCATED(MATRIX%DATA_INTG)) DEALLOCATE(MATRIX%DATA_INTG)
      IF(ALLOCATED(MATRIX%DATA_SP)) DEALLOCATE(MATRIX%DATA_SP)
      IF(ALLOCATED(MATRIX%DATA_DP)) DEALLOCATE(MATRIX%DATA_DP)
      IF(ALLOCATED(MATRIX%DATA_L)) DEALLOCATE(MATRIX%DATA_L)
      DEALLOCATE(MATRIX)
    ENDIF

    CALL EXITS("MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a matrix.
  SUBROUTINE MATRIX_INITIALISE(MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MATRIX_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      !!TODO: have a matrix user number etc.
      MATRIX%ID=0
      MATRIX%MATRIX_FINISHED=.FALSE.
      MATRIX%M=0
      MATRIX%N=0
      MATRIX%MAX_M=-1
      MATRIX%MAX_N=-1
      MATRIX%DATA_TYPE=0
      MATRIX%STORAGE_TYPE=0
      MATRIX%NUMBER_NON_ZEROS=0
      MATRIX%SIZE=0      
      MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW=0      
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("MATRIX_INITIALISE")
    RETURN 1
  END SUBROUTINE MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the maximum number of columns in each row of a distributed matrix.
  SUBROUTINE MATRIX_MAX_COLUMNS_PER_ROW_GET(MATRIX,MAX_COLUMNS_PER_ROW,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: MAX_COLUMNS_PER_ROW !<On return, the maximum number of columns in each row
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("MATRIX_MAX_COLUMNS_PER_ROW_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        MAX_COLUMNS_PER_ROW=MATRIX%MAXIMUM_COLUMN_INDICES_PER_ROW
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_MAX_COLUMNS_PER_ROW_GET")
    RETURN
999 CALL ERRORS("MATRIX_MAX_COLUMNS_PER_ROW_GET",ERR,ERROR)
    CALL EXITS("MATRIX_MAXCOLUMNS_PER_ROW_GET")
    RETURN 1
  END SUBROUTINE MATRIX_MAX_COLUMNS_PER_ROW_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of non zeros for a matrix.
  SUBROUTINE MATRIX_NUMBER_NON_ZEROS_SET(MATRIX,NUMBER_NON_ZEROS,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: NUMBER_NON_ZEROS !<The number of non zeros in the matrix to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_NUMBER_NON_ZEROS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("The matrix has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set the number of non-zeros for a matrix with block storage.",ERR,ERROR,*999)
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set the number of non-zeros for a matrix with diagonal storage.",ERR,ERROR,*999)
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set the number of non-zeros for a matrix with column major storage.",ERR,ERROR,*999)          
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set the number of non-zeros for a matrix with row major storage.",ERR,ERROR,*999)          
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
          IF(NUMBER_NON_ZEROS>=0) THEN
            MATRIX%NUMBER_NON_ZEROS=NUMBER_NON_ZEROS
          ELSE
            LOCAL_ERROR="The number of non-zeros ("//TRIM(NUMBER_TO_VSTRING(NUMBER_NON_ZEROS,"*",ERR,ERROR))// &
              & ") is invalid. The number must be greater than or equal to zero."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_NUMBER_NON_ZEROS_SET")
    RETURN
999 CALL ERRORS("MATRIX_NUMBER_NON_ZEROS_SET",ERR,ERROR)
    CALL EXITS("MATRIX_NUMBER_NON_ZEROS_SET")
    RETURN 1
  END SUBROUTINE MATRIX_NUMBER_NON_ZEROS_SET

  !
  !================================================================================================================================
  !

  !>Gets the number of non zeros for a matrix.
  SUBROUTINE MATRIX_NUMBER_NON_ZEROS_GET(MATRIX,NUMBER_NON_ZEROS,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: NUMBER_NON_ZEROS !<On return, The number of non zeros in the matrix to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_NUMBER_NON_ZEROS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE, &
          & MATRIX_ROW_MAJOR_STORAGE_TYPE,MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE, &
          & MATRIX_ROW_COLUMN_STORAGE_TYPE)
          NUMBER_NON_ZEROS=MATRIX%NUMBER_NON_ZEROS
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("The matrix is not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_NUMBER_NON_ZEROS_GET")
    RETURN
999 CALL ERRORS("MATRIX_NUMBER_NON_ZEROS_GET",ERR,ERROR)
    CALL EXITS("MATRIX_NUMBER_NON_ZEROS_GET")
    RETURN 1
  END SUBROUTINE MATRIX_NUMBER_NON_ZEROS_GET

  !
  !================================================================================================================================
  !
  ! TODO
  !>Sets the list of a matrix.
  SUBROUTINE MATRIX_LINKLIST_SET(MATRIX,LIST,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    TYPE(LinkedList),pointer :: LIST(:) 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("MATRIX_LINKLIST_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("The matrix has been finished",ERR,ERROR,*999)
      ELSE
        MATRIX%LIST => LIST      
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_LINKLIST_SET")
    RETURN
999 CALL ERRORS("MATRIX_LINKLIST_SET",ERR,ERROR)
    CALL EXITS("MATRIX_LINKLIST_SET")
    RETURN 1
  END SUBROUTINE MATRIX_LINKLIST_SET

  !
  !================================================================================================================================
  !
  ! TODO
  !>Gets the maximum number of columns in each row of a distributed matrix.
  SUBROUTINE MATRIX_LINKLIST_GET(MATRIX,LIST,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    type(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("MATRIX_LINKLIST_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        LIST=>MATRIX%LIST
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_LINKLIST_GET")
    RETURN
999 CALL ERRORS("MATRIX_LINKLIST_GET",ERR,ERROR)
    CALL EXITS("MATRIX_LINKLIST_GET")
    RETURN 1
  END SUBROUTINE MATRIX_LINKLIST_GET

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the maximum size of a matrix.
  SUBROUTINE MATRIX_MAX_SIZE_SET(MATRIX,MAX_M,MAX_N,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: MAX_M !<The maximum number of rows to set
    INTEGER(INTG), INTENT(IN) :: MAX_N !<The maximum number of columns to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_MAX_SIZE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("The matrix has been finished.",ERR,ERROR,*999)
      ELSE
        IF(MAX_M>0) THEN
          IF(MAX_N>0) THEN
            IF(MAX_M>=MATRIX%M) THEN
              IF(MAX_N>=MATRIX%N) THEN
                MATRIX%MAX_M=MAX_M
                MATRIX%MAX_N=MAX_N
              ELSE
                LOCAL_ERROR="The maximum number of matrix rows ("//TRIM(NUMBER_TO_VSTRING(MAX_N,"*",ERR,ERROR))// &
                  & ") must be >= the number of matrix rows ("//TRIM(NUMBER_TO_VSTRING(MATRIX%N,"*",ERR,ERROR))//")."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The maximum number of matrix columns ("//TRIM(NUMBER_TO_VSTRING(MAX_M,"*",ERR,ERROR))// &
                & ") must be >= the number of matrix columns ("//TRIM(NUMBER_TO_VSTRING(MATRIX%M,"*",ERR,ERROR))//")."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The maximum number of matrix columns of "//TRIM(NUMBER_TO_VSTRING(MAX_N,"*",ERR,ERROR))// &
              & " is invalid. The number must be > 0."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The maximum number of matrix rows of "//TRIM(NUMBER_TO_VSTRING(MAX_M,"*",ERR,ERROR))// &
            & " is invalid. The number must be > 0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_MAX_SIZE_SET")
    RETURN
999 CALL ERRORS("MATRIX_MAX_SIZE_SET",ERR,ERROR)
    CALL EXITS("MATRIX_MAX_SIZE_SET")
    RETURN 1
  END SUBROUTINE MATRIX_MAX_SIZE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the size of a matrix.
  SUBROUTINE MATRIX_OUTPUT(ID,MATRIX,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID to output to
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j
    CHARACTER(LEN=9) :: ROW_STRING,COL_STRING
    CHARACTER(LEN=39) :: INITIAL_STRING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          SELECT CASE(MATRIX%DATA_TYPE)
          CASE(MATRIX_VECTOR_INTG_TYPE)
            CALL WRITE_STRING_MATRIX(ID,1,1,MATRIX%M,1,1,MATRIX%N,8,8,RESHAPE(MATRIX%DATA_INTG,(/MATRIX%MAX_M,MATRIX%MAX_N/)), &
              & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,I13))','(20X,8(X,I13))', &
              & ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_SP_TYPE)
            CALL WRITE_STRING_MATRIX(ID,1,1,MATRIX%M,1,1,MATRIX%N,8,8,RESHAPE(MATRIX%DATA_SP,(/MATRIX%MAX_M,MATRIX%MAX_N/)), &
              & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,E13.6))','(20X,8(X,E13.6))', &
              & ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_DP_TYPE)
            CALL WRITE_STRING_MATRIX(ID,1,1,MATRIX%M,1,1,MATRIX%N,8,8,RESHAPE(MATRIX%DATA_DP,(/MATRIX%MAX_M,MATRIX%MAX_N/)), &
              & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,E13.6))','(20X,8(X,E13.6))', &
              & ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_L_TYPE)            
            CALL WRITE_STRING_MATRIX(ID,1,1,MATRIX%M,1,1,MATRIX%N,8,8,RESHAPE(MATRIX%DATA_L,(/MATRIX%MAX_M,MATRIX%MAX_N/)), &
              & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,L13))','(20X,8(X,L13))', &
              & ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The matrix data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          DO i=1,MATRIX%M
            ROW_STRING=NUMBER_TO_CHARACTER(i,"I9",ERR,ERROR)
            SELECT CASE(MATRIX%DATA_TYPE)              
            CASE(MATRIX_VECTOR_INTG_TYPE)
              INITIAL_STRING='("Matrix('//ROW_STRING//',:):",8(X,I13))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%ROW_INDICES(i),1,MATRIX%ROW_INDICES(i+1)-1,8,8,MATRIX%DATA_INTG,INITIAL_STRING, &
                & '(20X,8(X,I13))',ERR,ERROR,*999)
            CASE(MATRIX_VECTOR_SP_TYPE)
              INITIAL_STRING='("Matrix('//ROW_STRING//',:):",8(X,E13.6))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%ROW_INDICES(i),1,MATRIX%ROW_INDICES(i+1)-1,8,8,MATRIX%DATA_SP,INITIAL_STRING, &
                & '(20X,8(X,E13.6))',ERR,ERROR,*999)
            CASE(MATRIX_VECTOR_DP_TYPE)
              INITIAL_STRING='("Matrix('//ROW_STRING//',:):",8(X,E13.6))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%ROW_INDICES(i),1,MATRIX%ROW_INDICES(i+1)-1,8,8,MATRIX%DATA_DP,INITIAL_STRING, &
                & '(20X,8(X,E13.6))',ERR,ERROR,*999)
            CASE(MATRIX_VECTOR_L_TYPE)            
              INITIAL_STRING='("Matrix('//ROW_STRING//',:):",8(X,L13))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%ROW_INDICES(i),1,MATRIX%ROW_INDICES(i+1)-1,8,8,MATRIX%DATA_L,INITIAL_STRING, &
                & '(20X,8(X,L13))',ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The matrix data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDDO !i
        CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          DO j=1,MATRIX%N
            COL_STRING=NUMBER_TO_CHARACTER(j,"I9",ERR,ERROR)
            SELECT CASE(MATRIX%DATA_TYPE)              
            CASE(MATRIX_VECTOR_INTG_TYPE)
              INITIAL_STRING='("Matrix(:,'//COL_STRING//'):",8(X,I13))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%COLUMN_INDICES(j),1,MATRIX%COLUMN_INDICES(j+1)-1,8,8,MATRIX%DATA_INTG, &
                & INITIAL_STRING,'(20X,8(X,I13))',ERR,ERROR,*999)
            CASE(MATRIX_VECTOR_SP_TYPE)
              INITIAL_STRING='("Matrix(:,'//COL_STRING//'):",8(X,E13.6))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%COLUMN_INDICES(j),1,MATRIX%COLUMN_INDICES(j+1)-1,8,8,MATRIX%DATA_SP, &
                & INITIAL_STRING,'(20X,8(X,E13.6))',ERR,ERROR,*999)
            CASE(MATRIX_VECTOR_DP_TYPE)
              INITIAL_STRING='("Matrix(:,'//COL_STRING//'):",8(X,E13.6))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%COLUMN_INDICES(j),1,MATRIX%COLUMN_INDICES(j+1)-1,8,8,MATRIX%DATA_DP, &
                & INITIAL_STRING,'(20X,8(X,E13.6))',ERR,ERROR,*999)
            CASE(MATRIX_VECTOR_L_TYPE)            
              INITIAL_STRING='("Matrix(:,'//COL_STRING//'):",8(X,L13))'
              CALL WRITE_STRING_VECTOR(ID,MATRIX%COLUMN_INDICES(j),1,MATRIX%COLUMN_INDICES(j+1)-1,8,8,MATRIX%DATA_L, &
                & INITIAL_STRING,'(20X,8(X,L13))',ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The matrix data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDDO !j
        CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_OUTPUT")
    RETURN
999 CALL ERRORS("MATRIX_OUTPUT",ERR,ERROR)
    CALL EXITS("MATRIX_OUTPUT")
    RETURN 1
  END SUBROUTINE MATRIX_OUTPUT

  !
  !================================================================================================================================
  !

  !>Sets/changes the size of a matrix.
  SUBROUTINE MATRIX_SIZE_SET(MATRIX,M,N,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: M !<The number of rows to set
    INTEGER(INTG), INTENT(IN) :: N !<The number of columns to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_SIZE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("The matrix has been finished.",ERR,ERROR,*999)
      ELSE
        IF(M>0) THEN
          IF(N>0) THEN
            MATRIX%M=M
            MATRIX%N=N
          ELSE
            LOCAL_ERROR="The number of matrix columns of "//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))// &
              & " is invalid. The number must be >0."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The number of matrix rows of "//TRIM(NUMBER_TO_VSTRING(M,"*",ERR,ERROR))// &
            & " is invalid. The number must be >0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_SIZE_SET")
    RETURN
999 CALL ERRORS("MATRIX_SIZE_SET",ERR,ERROR)
    CALL EXITS("MATRIX_SIZE_SET")
    RETURN 1
  END SUBROUTINE MATRIX_SIZE_SET

  !
  !================================================================================================================================
  !

  !>Returns the storage location in the data array of a matrix that correponds to location I,J. If the location does not exist the routine returns zero.
  SUBROUTINE MATRIX_STORAGE_LOCATION_FIND(MATRIX,I,J,LOCATION,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: I !<The row number of the location to find
    INTEGER(INTG), INTENT(IN) :: J !<The column number of the location to find
    INTEGER(INTG), INTENT(OUT) :: LOCATION !<On return the location of the specified row and column in the matrix data. If the row and column does not exist in the matrix then zero is returned.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOWLIMIT,MIDPOINT,UPLIMIT
    LOGICAL :: FOUNDCOLUMN, FOUNDROW
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_STORAGE_LOCATION_FIND",ERR,ERROR,*999)

    LOCATION=0
    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(I<1.OR.I>MATRIX%M) THEN
          LOCAL_ERROR="Row number "//TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//" is outside the matrix range of 1 to "// &
            & TRIM(NUMBER_TO_VSTRING(MATRIX%M,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
        IF(J<1.OR.J>MATRIX%N) THEN
          LOCAL_ERROR="Column number "//TRIM(NUMBER_TO_VSTRING(J,"*",ERR,ERROR))//" is outside the matrix range of 1 to "// &
            & TRIM(NUMBER_TO_VSTRING(MATRIX%M,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          LOCATION=I+(J-1)*MATRIX%M
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          IF(I==J) LOCATION=I
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          LOCATION=I+(J-1)*MATRIX%MAX_M
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          LOCATION=(I-1)*MATRIX%MAX_N+J
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          !Search for the column number in the sparsity list using the bisection (binary search) algorithm
          LOWLIMIT=MATRIX%ROW_INDICES(I)
          IF(J>=MATRIX%COLUMN_INDICES(LOWLIMIT)) THEN
            UPLIMIT=MATRIX%ROW_INDICES(I+1)
            IF(UPLIMIT>LOWLIMIT) THEN
              IF(J<=MATRIX%COLUMN_INDICES(UPLIMIT-1)) THEN
                DO WHILE((UPLIMIT-LOWLIMIT)>1)
                  MIDPOINT=(UPLIMIT+LOWLIMIT)/2
                  IF(MATRIX%COLUMN_INDICES(MIDPOINT)>J) THEN
                    UPLIMIT=MIDPOINT
                  ELSE
                    LOWLIMIT=MIDPOINT
                  ENDIF
                ENDDO
                IF(MATRIX%COLUMN_INDICES(LOWLIMIT)==J) LOCATION=LOWLIMIT
              ENDIF
            ENDIF
          ENDIF
        CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          !Search for the row number in the sparsity list using the bisection (binary search) algorithm
          LOWLIMIT=MATRIX%COLUMN_INDICES(J)
          IF(I>=MATRIX%ROW_INDICES(LOWLIMIT)) THEN
            UPLIMIT=MATRIX%COLUMN_INDICES(J+1)
            IF(UPLIMIT>LOWLIMIT) THEN
              IF(I<=MATRIX%ROW_INDICES(UPLIMIT-1)) THEN
                DO WHILE((UPLIMIT-LOWLIMIT)>1)
                  MIDPOINT=(UPLIMIT+LOWLIMIT)/2
                  IF(MATRIX%ROW_INDICES(MIDPOINT)>I) THEN
                    UPLIMIT=MIDPOINT
                  ELSE
                    LOWLIMIT=MIDPOINT
                  ENDIF
                ENDDO
                IF(MATRIX%ROW_INDICES(LOWLIMIT)==I) LOCATION=LOWLIMIT
              ENDIF
            ENDIF
          ENDIF
        CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
          FOUNDROW=.FALSE.
          LOCATION=1
          DO WHILE(.NOT.FOUNDCOLUMN.AND.LOCATION<=MATRIX%SIZE)
            IF(MATRIX%ROW_INDICES(LOCATION)==I) THEN
              DO WHILE(.NOT.FOUNDCOLUMN.AND.LOCATION<=MATRIX%SIZE)
                IF(MATRIX%COLUMN_INDICES(LOCATION)==J.AND.MATRIX%ROW_INDICES(LOCATION)==I) THEN
                  FOUNDCOLUMN=.TRUE.
                ELSE IF(MATRIX%ROW_INDICES(LOCATION)/=I) THEN
                  LOCATION=MATRIX%SIZE+1
                ELSE
                  LOCATION=LOCATION+1
                ENDIF
              ENDDO
            ELSE
              LOCATION=LOCATION+1
            ENDIF
          ENDDO
          IF(.NOT.(FOUNDROW.AND.FOUNDCOLUMN)) LOCATION=0
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_STORAGE_LOCATION_FIND")
    RETURN
999 CALL ERRORS("MATRIX_STORAGE_LOCATION_FIND",ERR,ERROR)
    CALL EXITS("MATRIX_STORAGE_LOCATION_FIND")
    RETURN 1
  END SUBROUTINE MATRIX_STORAGE_LOCATION_FIND

  !
  !================================================================================================================================
  !

  !>Gets the storage locations (sparsity pattern) of a matrix.
  SUBROUTINE MATRIX_STORAGE_LOCATIONS_GET(MATRIX,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<ROW_INDICES(i). On return, the row index values for the matrix.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). On return, the column index values for the matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_STORAGE_LOCATIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not get matrix locations for a block storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not get matrix locations for a diagonal storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not get matrix locations for a column major storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not get matrix locations for a row major storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)          
          ROW_INDICES=>MATRIX%ROW_INDICES
          COLUMN_INDICES=>MATRIX%COLUMN_INDICES
        CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          ROW_INDICES=>MATRIX%ROW_INDICES
          COLUMN_INDICES=>MATRIX%COLUMN_INDICES          
        CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
          ROW_INDICES=>MATRIX%ROW_INDICES
          COLUMN_INDICES=>MATRIX%COLUMN_INDICES
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
   ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_STORAGE_LOCATIONS_GET")
    RETURN
999 CALL ERRORS("MATRIX_STORAGE_LOCATIONS_GET",ERR,ERROR)
    CALL EXITS("MATRIX_STORAGE_LOCATIONS_GET")
    RETURN 1
  END SUBROUTINE MATRIX_STORAGE_LOCATIONS_GET

  !
  !================================================================================================================================
  !

  !>Sets the storage locations (sparsity pattern) in a matrix to that specified by the row and column indices.
  SUBROUTINE MATRIX_STORAGE_LOCATIONS_SET(MATRIX,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index values for the sparisty pattern.
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column index values for the sparsity pattern.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_STORAGE_LOCATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("Matrix has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(MATRIX%STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set matrix locations for a block storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set matrix locations for a diagonal storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set matrix locations for a column major storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FLAG_ERROR("Can not set matrix locations for a row major storage matrix.",ERR,ERROR,*999)
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          IF(SIZE(ROW_INDICES,1)==MATRIX%M+1) THEN
            IF(SIZE(COLUMN_INDICES,1)==MATRIX%NUMBER_NON_ZEROS) THEN
              IF(ROW_INDICES(1)==1) THEN
                IF(ROW_INDICES(MATRIX%M+1)==MATRIX%NUMBER_NON_ZEROS+1) THEN
                  DO i=2,MATRIX%M+1
                    IF(ROW_INDICES(i)<ROW_INDICES(i-1)) THEN
                      LOCAL_ERROR="Invalid row indices. Row "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" index number ("// &
                        & TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//") is less than row "// &
                        & TRIM(NUMBER_TO_VSTRING(i-1,"*",ERR,ERROR))//" index number ("// &
                        & TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i-1),"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF                    
                  ENDDO !i
                  DO i=1,MATRIX%M
                    DO j=ROW_INDICES(i),ROW_INDICES(i+1)-1
                      k=COLUMN_INDICES(j)
                      IF(k>0) THEN
                        IF(k>MATRIX%N) THEN
                          LOCAL_ERROR="Invalid column indices. Column index "//TRIM(NUMBER_TO_VSTRING(j,"*",ERR,ERROR))//" ("// &
                            & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//") is greater than the number of columns ("// &
                            & TRIM(NUMBER_TO_VSTRING(MATRIX%N,"*",ERR,ERROR))//")."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Invalid column indices. Column index "//TRIM(NUMBER_TO_VSTRING(j,"*",ERR,ERROR))//" ("// &
                          & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//") is less than zero."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !j
                  ENDDO !i
                  IF(ALLOCATED(MATRIX%ROW_INDICES)) DEALLOCATE(MATRIX%ROW_INDICES)
                  IF(ALLOCATED(MATRIX%COLUMN_INDICES)) DEALLOCATE(MATRIX%COLUMN_INDICES)
                  ALLOCATE(MATRIX%ROW_INDICES(MATRIX%M+1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix row indices.",ERR,ERROR,*999)
                  ALLOCATE(MATRIX%COLUMN_INDICES(MATRIX%NUMBER_NON_ZEROS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix column indices.",ERR,ERROR,*999)                  
                  MATRIX%ROW_INDICES=ROW_INDICES
                  MATRIX%COLUMN_INDICES=COLUMN_INDICES
                  !Don't really need this???
                  !DO i=1,MATRIX%M
                  !  CALL LIST_SORT(MATRIX%COLUMN_INDICES(MATRIX%ROW_INDICES(i):MATRIX%ROW_INDICES(i+1)-1),ERR,ERROR,*999)
                  !ENDDO !i
                ELSE
                  LOCAL_ERROR="Invalid row indices. The last row index ("// &
                    & TRIM(NUMBER_TO_VSTRING(ROW_INDICES(MATRIX%M+1),"*",ERR,ERROR))// &
                    & ") does not equal the number of non-zeros + 1 ("// &
                    & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS+1,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Invalid row indices. The first row index ("// &
                  & TRIM(NUMBER_TO_VSTRING(ROW_INDICES(1),"*",ERR,ERROR))//") does not equal 1."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The supplied number of column indices ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
                & ") does not match the number of non-zeros in the matrix ("// &
                & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS,"*",ERR,ERROR))//")."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The supplied number of row indices ("//TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
              & ") does not match the number of rows in the matrix + 1 ("// &
              & TRIM(NUMBER_TO_VSTRING(MATRIX%M+1,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          IF(SIZE(COLUMN_INDICES,1)==MATRIX%N+1) THEN
            IF(SIZE(ROW_INDICES,1)==MATRIX%NUMBER_NON_ZEROS) THEN
             IF(COLUMN_INDICES(1)==1) THEN
                IF(COLUMN_INDICES(MATRIX%N+1)==MATRIX%NUMBER_NON_ZEROS+1) THEN
                  IF(COLUMN_INDICES(1)/=1) THEN
                    LOCAL_ERROR="Invalid column indices. Column index 1 ("// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(1),"*",ERR,ERROR))//") "// &
                      & " should be equal to one."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END IF
                  DO j=2,MATRIX%N+1
                    IF(COLUMN_INDICES(j)<COLUMN_INDICES(j-1)) THEN
                      LOCAL_ERROR="Invalid column indices. Column "//TRIM(NUMBER_TO_VSTRING(j,"*",ERR,ERROR))// &
                        & " index number ("//TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//") is less than column "// &
                        & TRIM(NUMBER_TO_VSTRING(j-1,"*",ERR,ERROR))//" index number ("// &
                        & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j-1),"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END IF
                    IF(COLUMN_INDICES(j)<0.OR.COLUMN_INDICES(j)>MATRIX%NUMBER_NON_ZEROS+1) THEN
                      LOCAL_ERROR="Invalid column indices. Column index "//TRIM(NUMBER_TO_VSTRING(j,"*",ERR,ERROR))//" ("// &
                        & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//") "// &
                        & " should be in the range of one to the number of non-zeros + 1 ("// &
                        & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS+1,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END IF
                  ENDDO !i
                  DO j=1,MATRIX%N
                    DO i=COLUMN_INDICES(j),COLUMN_INDICES(j+1)-1
                      k=ROW_INDICES(i)
                      IF(k>0) THEN
                        IF(k>MATRIX%M) THEN
                          LOCAL_ERROR="Invalid row indices. Row index "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" ("// &
                            & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//") is greater than the number of rows ("// &
                            & TRIM(NUMBER_TO_VSTRING(MATRIX%M,"*",ERR,ERROR))//")."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Invalid row indices. Row index "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" ("// &
                          & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//") is less than zero."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !i
                  ENDDO !j
                  IF(ALLOCATED(MATRIX%ROW_INDICES)) DEALLOCATE(MATRIX%ROW_INDICES)
                  IF(ALLOCATED(MATRIX%COLUMN_INDICES)) DEALLOCATE(MATRIX%COLUMN_INDICES)
                  ALLOCATE(MATRIX%ROW_INDICES(MATRIX%NUMBER_NON_ZEROS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix row indices.",ERR,ERROR,*999)
                  ALLOCATE(MATRIX%COLUMN_INDICES(MATRIX%N+1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix column indices.",ERR,ERROR,*999)                  
                  MATRIX%ROW_INDICES=ROW_INDICES
                  MATRIX%COLUMN_INDICES=COLUMN_INDICES
                  !Don't really need this???
                  !DO j=1,MATRIX%N                    
                  !  CALL LIST_SORT(MATRIX%ROW_INDICES(MATRIX%COLUMN_INDICES(j):MATRIX%COLUMN_INDICES(j+1)-1),ERR,ERROR,*999)
                  !ENDDO !j
                ELSE
                  LOCAL_ERROR="Invalid column indices. The last column index ("// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(MATRIX%N+1),"*",ERR,ERROR))// &
                    & ") does not equal the number of non-zeros + 1 ("// &
                    & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS+1,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Invalid column indices. The first column index ("// &
                  & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(1),"*",ERR,ERROR))//") does not equal 1."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The supplied number of row indices ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
                & ") does not match the number of non-zeros in the matrix ("// &
                & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS,"*",ERR,ERROR))//")."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The supplied number of column indices ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not match the number of columns in the matrix + 1 ("// &
              & TRIM(NUMBER_TO_VSTRING(MATRIX%N+1,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
          IF(SIZE(ROW_INDICES,1)==MATRIX%NUMBER_NON_ZEROS) THEN
            IF(SIZE(COLUMN_INDICES,1)==MATRIX%NUMBER_NON_ZEROS) THEN
              DO k=1,MATRIX%NUMBER_NON_ZEROS
                IF(ROW_INDICES(k)<1.OR.ROW_INDICES(k)>MATRIX%M) THEN
                  LOCAL_ERROR="Invalid row indices. Row index number "//TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" ("// &
                    & TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))// &
                    & ") is out of range. The row index must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(MATRIX%M,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE IF(COLUMN_INDICES(k)<1.OR.COLUMN_INDICES(k)>MATRIX%N) THEN
                  LOCAL_ERROR="Invalid column indices. Column index number "//TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" ("// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))// &
                    & ") is out of range. The column index must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(MATRIX%N,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !k
              MATRIX%ROW_INDICES=ROW_INDICES
              MATRIX%COLUMN_INDICES=COLUMN_INDICES
              !!TODO: sort the row and colum indices!!!!!
            ELSE
              LOCAL_ERROR="The supplied number of column indices ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
                & ") does not match the number of non-zeros in the matrix ("// &
                & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS,"*",ERR,ERROR))//")."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The supplied number of row indices ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
              & ") does not match the number of non-zeros in the matrix ("// &
              & TRIM(NUMBER_TO_VSTRING(MATRIX%NUMBER_NON_ZEROS,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_STORAGE_LOCATIONS_SET")
    RETURN
999 CALL ERRORS("MATRIX_STORAGE_LOCATIONS_SET",ERR,ERROR)
    CALL EXITS("MATRIX_STORAGE_LOCATIONS_SET")
    RETURN 1
  END SUBROUTINE MATRIX_STORAGE_LOCATIONS_SET

  !
  !================================================================================================================================
  !

  !>Gets the storage type for a matrix.
  SUBROUTINE MATRIX_STORAGE_TYPE_GET(MATRIX,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: STORAGE_TYPE !<On return, the storage type of the matrix. \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MATRIX_STORAGE_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN        
        STORAGE_TYPE=MATRIX%STORAGE_TYPE
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_STORAGE_TYPE_GET")
    RETURN
999 CALL ERRORS("MATRIX_STORAGE_TYPE_GET",ERR,ERROR)
    CALL EXITS("MATRIX_STORAGE_TYPE_GET")
    RETURN 1
  END SUBROUTINE MATRIX_STORAGE_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the storage type for a matrix.
  SUBROUTINE MATRIX_STORAGE_TYPE_SET(MATRIX,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE !<The storage type to set. \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        CALL FLAG_ERROR("The matrix has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(STORAGE_TYPE)
        CASE(MATRIX_BLOCK_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
        CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_DIAGONAL_STORAGE_TYPE
        CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
        CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_ROW_MAJOR_STORAGE_TYPE
        CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
        CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
        CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
          MATRIX%STORAGE_TYPE=MATRIX_ROW_COLUMN_STORAGE_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The matrix storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("MATRIX_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("MATRIX_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE MATRIX_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Adds values to an integer matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_INTG(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_INTG(LOCATION)=MATRIX%DATA_INTG(LOCATION)+VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_INTG")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_INTG",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_INTG")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_INTG

  !
  !================================================================================================================================
  !

  !>Adds a value to an integer matrix at the location specified by the row and column index i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_INTG1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index for the value to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_INTG(LOCATION)=MATRIX%DATA_INTG(LOCATION)+VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_INTG1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_INTG1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_INTG1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_INTG1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to an integer matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_INTG2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_INTG2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_INTG(LOCATION)=MATRIX%DATA_INTG(LOCATION)+VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_INTG2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_INTG2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_INTG2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_INTG2

  !
  !================================================================================================================================
  !

  !>Adds values to a single precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_SP(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to add
    REAL(SP), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_SP(LOCATION)=MATRIX%DATA_SP(LOCATION)+VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_SP")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_SP",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_SP")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_SP

  !
  !================================================================================================================================
  !

  !>Adds a value to a single precision real matrix at the location specified by the row and column index i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_SP1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index for the value to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_SP(LOCATION)=MATRIX%DATA_SP(LOCATION)+VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_SP1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_SP1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_SP1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_SP1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a single precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_SP2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to add
    REAL(SP), INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_SP2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_SP(LOCATION)=MATRIX%DATA_SP(LOCATION)+VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_SP2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_SP2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_SP2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_SP2

  !
  !================================================================================================================================
  !

  !>Adds values to a double precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_DP(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to add
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_DP(LOCATION)=MATRIX%DATA_DP(LOCATION)+VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_DP")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_DP",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_DP")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_DP

  !
  !================================================================================================================================
  !

  !>Adds a value to a double precision real matrix at the location specified by the row and column index i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_DP1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index for the value to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_DP(LOCATION)=MATRIX%DATA_DP(LOCATION)+VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_DP1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_DP1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_DP1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_DP1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a double precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J)+VALUE
  SUBROUTINE MATRIX_VALUES_ADD_DP2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to add
    REAL(DP), INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_DP2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_DP(LOCATION)=MATRIX%DATA_DP(LOCATION)+VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_DP2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_DP2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_DP2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_DP2

  !
  !================================================================================================================================
  !

  !>Adds values to a logical matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J).OR.VALUE
  SUBROUTINE MATRIX_VALUES_ADD_L(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to add
    LOGICAL, INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_L",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_L(LOCATION)=MATRIX%DATA_L(LOCATION).OR.VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_L")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_L",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_L")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_L

  !
  !================================================================================================================================
  !

  !>Adds a value to a logical matrix at the location specified by the row and column index i.e., MATRIX(I,J)=MATRIX(I,J).OR.VALUE
  SUBROUTINE MATRIX_VALUES_ADD_L1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index for the value to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_L1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_L(LOCATION)=MATRIX%DATA_L(LOCATION).OR.VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_L1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_L1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_L1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_L1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a logical matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=MATRIX(I,J).OR.VALUE
  SUBROUTINE MATRIX_VALUES_ADD_L2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to add
    LOGICAL, INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_ADD_L2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_L(LOCATION)=MATRIX%DATA_L(LOCATION).OR.VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_ADD_L2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_ADD_L2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_ADD_L2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_ADD_L2

  !
  !================================================================================================================================
  !

  !>Gets the values in an integer matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_INTG(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  VALUES(k)=0
                ELSE
                  VALUES(k)=MATRIX%DATA_INTG(LOCATION)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_INTG")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_INTG",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_INTG")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_INTG

  !
  !================================================================================================================================
  !

  !>Gets a value in an integer matrix at the location specified by the row and column index i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_INTG1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to get
    INTEGER(INTG), INTENT(OUT) :: VALUE !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            VALUE=0
          ELSE
            VALUE=MATRIX%DATA_INTG(LOCATION)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_INTG1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_INTG1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_INTG1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_INTG1

  !
  !================================================================================================================================
  !

  !>Gets the matrix of values in an integer matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_INTG2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: VALUES(:,:) !<VALUES(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_INTG2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    VALUES(i,j)=0
                  ELSE
                    VALUES(i,j)=MATRIX%DATA_INTG(LOCATION)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_VALUES_GET_INTG2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_INTG2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_INTG2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_INTG2

  !
  !================================================================================================================================
  !

  !>Gets values in a single precision real matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_SP(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to get
    REAL(SP), INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  VALUES(k)=0.0_SP
                ELSE
                  VALUES(k)=MATRIX%DATA_SP(LOCATION)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_SP")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_SP",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_SP")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_SP

  !
  !================================================================================================================================
  !

  !>Gets a value in a single precision real matrix at the location specified by the row and column index i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_SP1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to get
    REAL(SP), INTENT(OUT) :: VALUE !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            VALUE=0.0_SP
          ELSE
            VALUE=MATRIX%DATA_SP(LOCATION)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_SP1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_SP1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_SP1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_SP1

  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a single precision real matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_SP2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to get
    REAL(SP), INTENT(OUT) :: VALUES(:,:) !<VALUES(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_SP2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    VALUES(i,j)=0.0_SP
                  ELSE
                    VALUES(i,j)=MATRIX%DATA_SP(LOCATION)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_SP2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_SP2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_SP2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_SP2

  !
  !================================================================================================================================
  !

  !>Gets values in a double precision real matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_DP(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to get
    REAL(DP), INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  VALUES(k)=0.0_DP
                ELSE
                  VALUES(k)=MATRIX%DATA_DP(LOCATION)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_DP")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_DP",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_DP")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_DP

  !
  !================================================================================================================================
  !

  !>Gets a value in a double precision real matrix at the location specified by the row and column index i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_DP1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to get
    REAL(DP), INTENT(OUT) :: VALUE !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            VALUE=0.0_DP
          ELSE
            VALUE=MATRIX%DATA_DP(LOCATION)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_DP1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_DP1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_DP1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_DP1

  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a double precision real matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_DP2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to get
    REAL(DP), INTENT(OUT) :: VALUES(:,:) !<VALUES(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_DP2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    VALUES(i,j)=0.0_DP
                  ELSE
                    VALUES(i,j)=MATRIX%DATA_DP(LOCATION)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_DP2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_DP2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_DP2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_DP2

  !
  !================================================================================================================================
  !

  !>Gets values in a logical matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_L(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to get
    LOGICAL, INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  VALUES(k)=.FALSE.
                ELSE
                  VALUES(k)=MATRIX%DATA_L(LOCATION)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_L")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_L",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_L")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_L

  !
  !================================================================================================================================
  !

  !>Gets a value in a logical matrix at the location specified by the row and column index i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_L1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to get
    LOGICAL, INTENT(OUT) :: VALUE !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_L1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            VALUE=.FALSE.
          ELSE
            VALUE=MATRIX%DATA_L(LOCATION)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_L1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_L1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_L1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_L1

  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a logical matrix at the location specified by the row and column indices i.e., VALUE=MATRIX(I,J)
  SUBROUTINE MATRIX_VALUES_GET_L2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to get
    LOGICAL, INTENT(OUT) :: VALUES(:,:) !<VALUES(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_GET_L2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    VALUES(i,j)=.FALSE.
                  ELSE
                    VALUES(i,j)=MATRIX%DATA_L(LOCATION)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_GET_L2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_GET_L2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_GET_L2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_GET_L2

  !
  !================================================================================================================================
  !

  !>Sets the values in an integer matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_INTG(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(i). The column index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_INTG(LOCATION)=VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_INTG")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_INTG",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_INTG")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_INTG

  !
  !================================================================================================================================
  !

  !>Sets a value in an integer matrix at the location specified by the row and column index i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_INTG1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to set
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_INTG(LOCATION)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_INTG1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_INTG1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_INTG1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_INTG1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in an integer matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_INTG2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INIDICES(j). The column index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the i,j'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_INTG2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_INTG(LOCATION)=VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_VALUES_SET_INTG2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_INTG2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_INTG2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_INTG2
  
  !
  !================================================================================================================================
  !

  !>Sets the values in a single precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_SP(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to set
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index of the i'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column index of the i'th value to set
    REAL(SP), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_SP(LOCATION)=VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_SP")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_SP",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_SP")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_SP

  !
  !================================================================================================================================
  !

  !>Sets the value in a single precision real matrix at the location specified by the row and column index i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_SP1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to set
    REAL(SP), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_SP(LOCATION)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_SP1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_SP1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_SP1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_SP1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in a single precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_SP2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to set
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index of the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(j). The column index of the ij'th value to set
    REAL(SP), INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_SP2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_SP(LOCATION)=VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MATRIX_VALUES_SET_SP2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_SP2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_SP2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_SP2

  !
  !================================================================================================================================
  !

  !>Sets the values in a double precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_DP(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to set.
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index of the i'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column index of the i'th value to set
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_DP(LOCATION)=VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_DP")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_DP",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_DP")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_DP

  !
  !================================================================================================================================
  !

  !>Sets a value in a double precision real matrix at the location specified by the row and column index i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_DP1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to set
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_DP(LOCATION)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_DP1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_DP1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_DP1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_DP1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in a double precision real matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_DP2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to set.
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index of the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(j). The column index of the ij'th value to set
    REAL(DP), INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_DP2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(COLUMN_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_DP(LOCATION)=VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_DP2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_DP2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_DP2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_DP2

  !
  !================================================================================================================================
  !

  !>Sets the values in a logical matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_L(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to set
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index of the i'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(i). The column index of the i'th value to set
    LOGICAL, INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: k,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,1)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
              DO k=1,SIZE(ROW_INDICES,1)
                CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(k),COLUMN_INDICES(k),LOCATION,ERR,ERROR,*999)
                IF(LOCATION==0) THEN
                  LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(k),"*",ERR,ERROR))//" and column "// &
                    & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(k),"*",ERR,ERROR))//" does not exist in the matrix."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  MATRIX%DATA_L(LOCATION)=VALUES(k)
                ENDIF
              ENDDO !k
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_L")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_L",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_L")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_L

  !
  !================================================================================================================================
  !

  !>Sets a value in a logical matrix at the location specified by the row and column index i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_L1(MATRIX,ROW_INDEX,COLUMN_INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: ROW_INDEX !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDEX !<The column index of the value to set
    LOGICAL, INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_L1",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDEX,COLUMN_INDEX,LOCATION,ERR,ERROR,*999)
          IF(LOCATION==0) THEN
            LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDEX,"*",ERR,ERROR))//" and column "// &
              & TRIM(NUMBER_TO_VSTRING(COLUMN_INDEX,"*",ERR,ERROR))//" does not exist in the matrix."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            MATRIX%DATA_L(LOCATION)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_L1")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_L1",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_L1")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_L1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in a logical matrix at the location specified by the row and column indices i.e., MATRIX(I,J)=VALUE
  SUBROUTINE MATRIX_VALUES_SET_L2(MATRIX,ROW_INDICES,COLUMN_INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(MATRIX_TYPE), POINTER :: MATRIX !<A pointer to the matrix to set
    INTEGER(INTG), INTENT(IN) :: ROW_INDICES(:) !<ROW_INDICES(i). The row index of the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: COLUMN_INDICES(:) !<COLUMN_INDICES(j). The column index of the ij'th value to set
    LOGICAL, INTENT(IN) :: VALUES(:,:) !<VALUES(i,j). The value of the ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LOCATION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MATRIX_VALUES_SET_L2",ERR,ERROR,*999)

    IF(ASSOCIATED(MATRIX)) THEN
      IF(MATRIX%MATRIX_FINISHED) THEN
        IF(SIZE(ROW_INDICES,1)==SIZE(VALUES,1)) THEN
          IF(SIZE(COLUMN_INDICES,1)==SIZE(VALUES,2)) THEN
            IF(MATRIX%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
              DO i=1,SIZE(ROW_INDICES,1)
                DO j=1,SIZE(ROW_INDICES,1)
                  CALL MATRIX_STORAGE_LOCATION_FIND(MATRIX,ROW_INDICES(i),COLUMN_INDICES(j),LOCATION,ERR,ERROR,*999)
                  IF(LOCATION==0) THEN
                    LOCAL_ERROR="Row "//TRIM(NUMBER_TO_VSTRING(ROW_INDICES(i),"*",ERR,ERROR))//" and column "// &
                      & TRIM(NUMBER_TO_VSTRING(COLUMN_INDICES(j),"*",ERR,ERROR))//" does not exist in the matrix."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ELSE
                    MATRIX%DATA_L(LOCATION)=VALUES(i,j)
                  ENDIF
                ENDDO !j
              ENDDO !i
            ELSE
              LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(MATRIX%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given values."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The size of the column indices array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(COLUMN_INDICES,1),"*",ERR,ERROR))// &
              & ") does not conform to the number of columns in the values array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,2),"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the row indices array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the number of rows in the values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The matrix has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MATRIX_VALUES_SET_L2")
    RETURN
999 CALL ERRORS("MATRIX_VALUES_SET_L2",ERR,ERROR)
    CALL EXITS("MATRIX_VALUES_SET_L2")
    RETURN 1
  END SUBROUTINE MATRIX_VALUES_SET_L2

  !
  !================================================================================================================================
  !

  !>Sets all values in an integer vector to the specified value.
  SUBROUTINE VECTOR_ALL_VALUES_SET_INTG(VECTOR,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector to set
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_ALL_VALUES_SET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          VECTOR%DATA_INTG=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_ALL_VALUES_SET_INTG")
    RETURN
999 CALL ERRORS("VECTOR_ALL_VALUES_SET_INTG",ERR,ERROR)
    CALL EXITS("VECTOR_ALL_VALUES_SET_INTG")
    RETURN 1
  END SUBROUTINE VECTOR_ALL_VALUES_SET_INTG

  !
  !================================================================================================================================
  !

  !>Sets all values in a single precision vector to the specified value.
  SUBROUTINE VECTOR_ALL_VALUES_SET_SP(VECTOR,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector to set
    REAL(SP), INTENT(IN) :: VALUE !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_ALL_VALUES_SET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          VECTOR%DATA_SP=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_ALL_VALUES_SET_SP")
    RETURN
999 CALL ERRORS("VECTOR_ALL_VALUES_SET_SP",ERR,ERROR)
    CALL EXITS("VECTOR_ALL_VALUES_SET_SP")
    RETURN 1
  END SUBROUTINE VECTOR_ALL_VALUES_SET_SP

  !
  !================================================================================================================================
  !

  !>Sets all values in a double precision vector to the specified value.
  SUBROUTINE VECTOR_ALL_VALUES_SET_DP(VECTOR,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector to set
    REAL(DP), INTENT(IN) :: VALUE !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_ALL_VALUES_SET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          VECTOR%DATA_DP=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_ALL_VALUES_SET_DP")
    RETURN
999 CALL ERRORS("VECTOR_ALL_VALUES_SET_DP",ERR,ERROR)
    CALL EXITS("VECTOR_ALL_VALUES_SET_DP")
    RETURN 1
  END SUBROUTINE VECTOR_ALL_VALUES_SET_DP

  !
  !================================================================================================================================
  !

  !> Sets all values in a logical vector to the specified value.
  SUBROUTINE VECTOR_ALL_VALUES_SET_L(VECTOR,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector to set
    LOGICAL, INTENT(IN) :: VALUE !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_ALL_VALUES_SET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          VECTOR%DATA_L=VALUE
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_ALL_VALUES_SET_L")
    RETURN
999 CALL ERRORS("VECTOR_ALL_VALUES_SET_L",ERR,ERROR)
    CALL EXITS("VECTOR_ALL_VALUES_SET_L")
    RETURN 1
  END SUBROUTINE VECTOR_ALL_VALUES_SET_L

  !
  !================================================================================================================================
  !

  !>Finihses the creation of a vector. 
  SUBROUTINE VECTOR_CREATE_FINISH(VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        CALL FLAG_ERROR("Vector has been finished.",ERR,ERROR,*999)
      ELSE
        IF(VECTOR%SIZE>0) THEN
          SELECT CASE(VECTOR%DATA_TYPE)
          CASE(MATRIX_VECTOR_INTG_TYPE)
            ALLOCATE(VECTOR%DATA_INTG(VECTOR%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate vector integer data.",ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_SP_TYPE)
            ALLOCATE(VECTOR%DATA_SP(VECTOR%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate vector single precision data.",ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_DP_TYPE)
            ALLOCATE(VECTOR%DATA_DP(VECTOR%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate vector double precision data.",ERR,ERROR,*999)
          CASE(MATRIX_VECTOR_L_TYPE)
            ALLOCATE(VECTOR%DATA_L(VECTOR%SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate vector logical data.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The vector data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
        VECTOR%ID=MATRIX_VECTOR_ID
        MATRIX_VECTOR_ID=MATRIX_VECTOR_ID+1
        VECTOR%VECTOR_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("VECTOR_CREATE_FINISH")
    RETURN
999 CALL ERRORS("VECTOR_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("VECTOR_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE VECTOR_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation a vector. 
  SUBROUTINE VECTOR_CREATE_START(VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("VECTOR_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      CALL FLAG_ERROR("Vector is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(VECTOR,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate the vector.",ERR,ERROR,*999)
      CALL VECTOR_INITIALISE(VECTOR,ERR,ERROR,*999)
      !Set the defaults
      VECTOR%DATA_TYPE=MATRIX_VECTOR_DP_TYPE
    ENDIF
 
    CALL EXITS("VECTOR_CREATE_START")
    RETURN
999 IF(ASSOCIATED(VECTOR)) CALL VECTOR_FINALISE(VECTOR,ERR,ERROR,*998)
998 CALL ERRORS("VECTOR_CREATE_START",ERR,ERROR)
    CALL EXITS("VECTOR_CREATE_START")
    RETURN 1
  END SUBROUTINE VECTOR_CREATE_START

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of an integer vector. Note: the values can be used for read operations but a VECTOR_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE VECTOR_DATA_GET_INTG(VECTOR,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), POINTER :: DATA(:) !<On return a pointer to the vector data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_DATA_GET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(VECTOR%VECTOR_FINISHED) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
            DATA=>VECTOR%DATA_INTG
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the integer data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("VECTOR_DATA_GET_INTG")
    RETURN
999 CALL ERRORS("VECTOR_DATA_GET_INTG",ERR,ERROR)
    CALL EXITS("VECTOR_DATA_GET_INTG")
    RETURN 1
  END SUBROUTINE VECTOR_DATA_GET_INTG

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a single precision vector. Note: the values can be used for read operations but a VECTOR_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE VECTOR_DATA_GET_SP(VECTOR,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    REAL(SP), POINTER :: DATA(:) !<On return a pointer to the vector data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_DATA_GET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(VECTOR%VECTOR_FINISHED) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
            DATA=>VECTOR%DATA_SP
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the single precision data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("VECTOR_DATA_GET_SP")
    RETURN
999 CALL ERRORS("VECTOR_DATA_GET_SP",ERR,ERROR)
    CALL EXITS("VECTOR_DATA_GET_SP")
    RETURN 1
  END SUBROUTINE VECTOR_DATA_GET_SP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a double precision vector. Note: the values can be used for read operations but a VECTOR_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE VECTOR_DATA_GET_DP(VECTOR,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    REAL(DP), POINTER :: DATA(:) !<On return a pointer to the vector data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_DATA_GET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(VECTOR%VECTOR_FINISHED) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
            DATA=>VECTOR%DATA_DP
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the double precision data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("VECTOR_DATA_GET_DP")
    RETURN
999 CALL ERRORS("VECTOR_DATA_GET_DP",ERR,ERROR)
    CALL EXITS("VECTOR_DATA_GET_DP")
    RETURN 1
  END SUBROUTINE VECTOR_DATA_GET_DP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a logical vector. Note: the values can be used for read operations but a VECTOR_VALUES_SET call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE VECTOR_DATA_GET_L(VECTOR,DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    LOGICAL, POINTER :: DATA(:) !<On return a pointer to the vector data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_DATA_GET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(ASSOCIATED(DATA)) THEN
        CALL FLAG_ERROR("Data is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DATA)
        IF(VECTOR%VECTOR_FINISHED) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
            DATA=>VECTOR%DATA_L
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the logical data type of the requested values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("VECTOR_DATA_GET_L")
    RETURN
999 CALL ERRORS("VECTOR_DATA_GET_L",ERR,ERROR)
    CALL EXITS("VECTOR_DATA_GET_L")
    RETURN 1
  END SUBROUTINE VECTOR_DATA_GET_L

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type of a vector.
  SUBROUTINE VECTOR_DATA_TYPE_SET(VECTOR,DATA_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector.
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type to set. \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_DATA_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        CALL FLAG_ERROR("The vector has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DATA_TYPE)
        CASE(MATRIX_VECTOR_INTG_TYPE)
          VECTOR%DATA_TYPE=MATRIX_VECTOR_INTG_TYPE
        CASE(MATRIX_VECTOR_SP_TYPE)
          VECTOR%DATA_TYPE=MATRIX_VECTOR_SP_TYPE
        CASE(MATRIX_VECTOR_DP_TYPE)
          VECTOR%DATA_TYPE=MATRIX_VECTOR_DP_TYPE
        CASE(MATRIX_VECTOR_L_TYPE)
          VECTOR%DATA_TYPE=MATRIX_VECTOR_L_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The vector data type of "//TRIM(NUMBER_TO_VSTRING(DATA_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_DATA_TYPE_SET")
    RETURN
999 CALL ERRORS("VECTOR_DATA_TYPE_SET",ERR,ERROR)
    CALL EXITS("VECTOR_DATA_TYPE_SET")
    RETURN 1
  END SUBROUTINE VECTOR_DATA_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Destroys a vector
  SUBROUTINE VECTOR_DESTROY(VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("VECTOR_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      CALL VECTOR_FINALISE(VECTOR,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_DESTROY")
    RETURN
999 CALL ERRORS("VECTOR_DESTROY",ERR,ERROR)
    CALL EXITS("VECTOR_DESTROY")
    RETURN 1
  END SUBROUTINE VECTOR_DESTROY

  !
  !================================================================================================================================
  !

  !>Duplicates a vector structure and returns a pointer to the new vector in NEW_VECTOR.
  SUBROUTINE VECTOR_DUPLICATE(VECTOR,NEW_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector to duplicate
    TYPE(VECTOR_TYPE), POINTER :: NEW_VECTOR !<On return a pointer to the new duplicated vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("VECTOR_DUPLICATE",ERR,ERROR,*998)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(ASSOCIATED(NEW_VECTOR)) THEN
        CALL FLAG_ERROR("New vector is already associated.",ERR,ERROR,*998)
      ELSE
        CALL VECTOR_CREATE_START(NEW_VECTOR,ERR,ERROR,*999)
        CALL VECTOR_DATA_TYPE_SET(NEW_VECTOR,VECTOR%DATA_TYPE,ERR,ERROR,*999)
        CALL VECTOR_SIZE_SET(NEW_VECTOR,VECTOR%N,ERR,ERROR,*999)
        CALL VECTOR_CREATE_FINISH(NEW_VECTOR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("VECTOR_DUPLICATE")
    RETURN
999 CALL VECTOR_FINALISE(NEW_VECTOR,ERR,ERROR,*998)
998 CALL ERRORS("VECTOR_DUPLICATE",ERR,ERROR)
    CALL EXITS("VECTOR_DUPLICATE")
    RETURN 1
  END SUBROUTINE VECTOR_DUPLICATE

  !
  !================================================================================================================================
  !

  !>Finalises a vector and deallocates all memory
  SUBROUTINE VECTOR_FINALISE(VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("VECTOR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(ALLOCATED(VECTOR%DATA_INTG)) DEALLOCATE(VECTOR%DATA_INTG)
      IF(ALLOCATED(VECTOR%DATA_SP)) DEALLOCATE(VECTOR%DATA_SP)
      IF(ALLOCATED(VECTOR%DATA_DP)) DEALLOCATE(VECTOR%DATA_DP)
      IF(ALLOCATED(VECTOR%DATA_L)) DEALLOCATE(VECTOR%DATA_L)
      DEALLOCATE(VECTOR)
    ENDIF

    CALL EXITS("VECTOR_FINALISE")
    RETURN
999 CALL ERRORS("VECTOR_FINALISE",ERR,ERROR)
    CALL EXITS("VECTOR_FINALISE")
    RETURN 1
  END SUBROUTINE VECTOR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a vector
  SUBROUTINE VECTOR_INITIALISE(VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("VECTOR_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      !!TODO: have a vector user number etc.
      VECTOR%ID=0
      VECTOR%VECTOR_FINISHED=.FALSE.
      VECTOR%N=0
      VECTOR%DATA_TYPE=0
      VECTOR%SIZE=0      
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_INITIALISE")
    RETURN
999 CALL ERRORS("VECTOR_INITIALISE",ERR,ERROR)
    CALL EXITS("VECTOR_INITIALISE")
    RETURN 1
  END SUBROUTINE VECTOR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the size of a vector
  SUBROUTINE VECTOR_SIZE_SET(VECTOR,N,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: N !<The size of the vector to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_SIZE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        CALL FLAG_ERROR("The vector has been finished.",ERR,ERROR,*999)
      ELSE
        IF(N>0) THEN
          VECTOR%N=N
        ELSE
          LOCAL_ERROR="The size of the vector ("//TRIM(NUMBER_TO_VSTRING(N,"*",ERR,ERROR))// &
            & ") is invalid. The number must be >0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_SIZE_SET")
    RETURN
999 CALL ERRORS("VECTOR_SIZE_SET",ERR,ERROR)
    CALL EXITS("VECTOR_SIZE_SET")
    RETURN 1
  END SUBROUTINE VECTOR_SIZE_SET

  !
  !================================================================================================================================
  !

  !>Gets the values in an integer vector at the indices specified.
  SUBROUTINE VECTOR_VALUES_GET_INTG(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to get
    INTEGER(INTG), INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VALUES(i)=VECTOR%DATA_INTG(k)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the integer data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_INTG")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_INTG",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_INTG")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_INTG

  !
  !================================================================================================================================
  !

  !>Gets a value in an integer vector at the location specified by the index
  SUBROUTINE VECTOR_VALUES_GET_INTG1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index of the vector to get
    INTEGER(INTG), INTENT(OUT) :: VALUE !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VALUE=VECTOR%DATA_INTG(INDEX)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_INTG1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_INTG1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_INTG1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_INTG1

  !
  !================================================================================================================================
  !

  !>Gets the values in a single precision real vector at the indices specified
  SUBROUTINE VECTOR_VALUES_GET_SP(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to get
    REAL(SP), INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VALUES(i)=VECTOR%DATA_SP(k)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the single precision data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_SP")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_SP",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_SP")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_SP

  !
  !================================================================================================================================
  !

  !>Gets a value in a single precision vector at the location specified by the index
  SUBROUTINE VECTOR_VALUES_GET_SP1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index of the vector to get
    REAL(SP), INTENT(OUT) :: VALUE !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VALUE=VECTOR%DATA_SP(INDEX)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_SP1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_SP1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_SP1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_SP1

  !
  !================================================================================================================================
  !

  !>Gets the values in a double precision real vector at the indices specified.
  SUBROUTINE VECTOR_VALUES_GET_DP(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to get
    REAL(DP), INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VALUES(i)=VECTOR%DATA_DP(k)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the double precision data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_DP")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_DP",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_DP")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_DP

  !
  !================================================================================================================================
  !

  !>Gets a value in a double precision vector at the location specified by the index
  SUBROUTINE VECTOR_VALUES_GET_DP1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index of the vector to get
    REAL(DP), INTENT(OUT) :: VALUE !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VALUE=VECTOR%DATA_DP(INDEX)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_DP1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_DP1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_DP1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_DP1

  !
  !================================================================================================================================
  !

  !>Gets the values in a logical real vector at the indices specified.
  SUBROUTINE VECTOR_VALUES_GET_L(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to get
    LOGICAL, INTENT(OUT) :: VALUES(:) !<VALUES(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VALUES(i)=VECTOR%DATA_L(k)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the logical data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_L")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_L",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_L")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_L

  !
  !================================================================================================================================
  !

  !>Gets a value in a logical vector at the location specified by the index
  SUBROUTINE VECTOR_VALUES_GET_L1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index of the vector to get
    LOGICAL, INTENT(OUT) :: VALUE !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_GET_L1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VALUE=VECTOR%DATA_L(INDEX)
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_GET_L1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_GET_L1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_GET_L1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_GET_L1

  !
  !================================================================================================================================
  !

  !>Sets the values in an integer vector at the specified indices.
  SUBROUTINE VECTOR_VALUES_SET_INTG(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to set
    INTEGER(INTG), INTENT(IN) :: VALUES(:) !<VALUES(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VECTOR%DATA_INTG(k)=VALUES(i)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the integer data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_INTG")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_INTG",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_INTG")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_INTG

  !
  !================================================================================================================================
  !
  
  !>Sets a value in an integer vector at the specified index.
  SUBROUTINE VECTOR_VALUES_SET_INTG1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index to set
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_INTG_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VECTOR%DATA_INTG(INDEX)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the integer data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_INTG1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_INTG1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_INTG1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_INTG1

  !
  !================================================================================================================================
  !

  !>Sets the values in a single precision vector at the specified indices.
  SUBROUTINE VECTOR_VALUES_SET_SP(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to set
    REAL(SP), INTENT(IN) :: VALUES(:) !<VALUES(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
               VECTOR%DATA_SP(k)=VALUES(i)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the single precision data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_SP")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_SP",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_SP")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_SP

  !
  !================================================================================================================================
  !

  !>Sets a value in a single precision vector at the specified index.
  SUBROUTINE VECTOR_VALUES_SET_SP1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index to set
    REAL(SP), INTENT(IN) :: VALUE !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_SP_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VECTOR%DATA_SP(INDEX)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the single precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_SP1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_SP1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_SP1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_SP1

  !
  !================================================================================================================================
  !

  !>Sets the values in a double precision vector at the specified indices.
  SUBROUTINE VECTOR_VALUES_SET_DP(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to set
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VECTOR%DATA_DP(k)=VALUES(i)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the double precision data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_DP")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_DP",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_DP")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_DP

  !
  !================================================================================================================================
  !

  !>Sets a value in a double precision vector at the specified index.
  SUBROUTINE VECTOR_VALUES_SET_DP1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index to set
    REAL(DP), INTENT(IN) :: VALUE !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_DP_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VECTOR%DATA_DP(INDEX)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the double precision data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_DP1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_DP1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_DP1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_DP1

  !
  !================================================================================================================================
  !

  !>Sets the values in a logical vector at the specified indices.
  SUBROUTINE VECTOR_VALUES_SET_L(VECTOR,INDICES,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDICES(:) !<INDICES(i). The i'th index to set
    LOGICAL, INTENT(IN) :: VALUES(:) !<VALUES(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(SIZE(INDICES,1)==SIZE(VALUES,1)) THEN
          IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
            DO i=1,SIZE(INDICES,1)
              k=INDICES(i)
              IF(k<1.OR.k>VECTOR%N) THEN
                LOCAL_ERROR="Index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is invalid. The index is "// &
                    & TRIM(NUMBER_TO_VSTRING(k,"*",ERR,ERROR))//" and it must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                VECTOR%DATA_L(k)=VALUES(i)
              ENDIF
            ENDDO !i
          ELSE
            LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not correspond to the logical data type of the given values."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The size of the indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(INDICES,1),"*",ERR,ERROR))// &
            & ") does not conform to the size of the values array ("//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_L")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_L",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_L")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_L

  !
  !================================================================================================================================
  !

  !>Sets a value in a logical vector at the specified index.
  SUBROUTINE VECTOR_VALUES_SET_L1(VECTOR,INDEX,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(VECTOR_TYPE), POINTER :: VECTOR !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: INDEX !<The index to set
    LOGICAL, INTENT(IN) :: VALUE !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("VECTOR_VALUES_SET_L1",ERR,ERROR,*999)

    IF(ASSOCIATED(VECTOR)) THEN
      IF(VECTOR%VECTOR_FINISHED) THEN
        IF(VECTOR%DATA_TYPE==MATRIX_VECTOR_L_TYPE) THEN
          IF(INDEX<1.OR.INDEX>VECTOR%N) THEN
            LOCAL_ERROR="The index value of "//TRIM(NUMBER_TO_VSTRING(INDEX,"*",ERR,ERROR))// &
              & " is invalid. The index must be between 1 and  "//TRIM(NUMBER_TO_VSTRING(VECTOR%N,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            VECTOR%DATA_L(INDEX)=VALUE
          ENDIF
        ELSE
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(VECTOR%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not correspond to the logical data type of the given value."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The vector has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Vector is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("VECTOR_VALUES_SET_L1")
    RETURN
999 CALL ERRORS("VECTOR_VALUES_SET_L1",ERR,ERROR)
    CALL EXITS("VECTOR_VALUES_SET_L1")
    RETURN 1
  END SUBROUTINE VECTOR_VALUES_SET_L1

  !
  !================================================================================================================================
  !
  
END MODULE MATRIX_VECTOR
