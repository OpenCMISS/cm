!> \file
!> \author Chris Bradley
!> \brief This module is a CMISS buffer module to the PETSc library.
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

!> This module is a CMISS buffer module to the PETSc library.
MODULE CMISS_PETSC
  
  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE KINDS
  USE ISO_VARYING_STRING
  USE TYPES
  
  IMPLICIT NONE
 
  PRIVATE

#include "petscversion.h"
  
#include "finclude/petsc.h"
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 1 )
#include "finclude/petscis.h"
#include "finclude/petscksp.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscsnes.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#endif
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
#include "finclude/petscts.h"
#endif

  !Module parameters

  !Insert mode types
  InsertMode, PARAMETER :: PETSC_INSERT_VALUES = INSERT_VALUES
  InsertMode, PARAMETER :: PETSC_ADD_VALUES = ADD_VALUES

  !Scatter mode types
  ScatterMode, PARAMETER :: PETSC_SCATTER_FORWARD = SCATTER_FORWARD
  ScatterMode, PARAMETER :: PETSC_SCATTER_REVERSE = SCATTER_REVERSE
  
  !KSPConvergedReason types
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_RTOL = KSP_CONVERGED_RTOL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_ATOL = KSP_CONVERGED_ATOL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_ITS = KSP_CONVERGED_ITS
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_CG_NEG_CURVE = KSP_CONVERGED_CG_NEG_CURVE
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_CG_CONSTRAINED = KSP_CONVERGED_CG_CONSTRAINED
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_STEP_LENGTH = KSP_CONVERGED_STEP_LENGTH
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_HAPPY_BREAKDOWN = KSP_CONVERGED_HAPPY_BREAKDOWN
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NULL = KSP_DIVERGED_NULL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_ITS = KSP_DIVERGED_ITS
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_DTOL = KSP_DIVERGED_DTOL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_BREAKDOWN = KSP_DIVERGED_BREAKDOWN
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_BREAKDOWN_BICG = KSP_DIVERGED_BREAKDOWN_BICG
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NONSYMMETRIC = KSP_DIVERGED_NONSYMMETRIC
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_INDEFINITE_PC = KSP_DIVERGED_INDEFINITE_PC
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 3 )
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NAN = KSP_DIVERGED_NAN
#else
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NAN = KSP_DIVERGED_NANORINF
#endif
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_INDEFINITE_MAT = KSP_DIVERGED_INDEFINITE_MAT
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_ITERATING = KSP_CONVERGED_ITERATING  
  
  !KSP types
  KSPType, PARAMETER :: PETSC_KSPRICHARDSON = KSPRICHARDSON
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  KSPType, PARAMETER :: PETSC_KSPCHEBYSHEV = KSPCHEBYSHEV
#else
  KSPType, PARAMETER :: PETSC_KSPCHEBYCHEV = KSPCHEBYCHEV
#endif  
  KSPType, PARAMETER :: PETSC_KSPCG = KSPCG
  KSPType, PARAMETER :: PETSC_KSPCGNE = KSPCGNE
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  KSPType, PARAMETER :: PETSC_KSPNASH = KSPNASH
#endif
  KSPType, PARAMETER :: PETSC_KSPSTCG = KSPSTCG
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  KSPType, PARAMETER :: PETSC_KSPGLTR = KSPGLTR
#endif
  KSPType, PARAMETER :: PETSC_KSPGMRES = KSPGMRES
  KSPType, PARAMETER :: PETSC_KSPFGMRES = KSPFGMRES
  KSPType, PARAMETER :: PETSC_KSPLGMRES = KSPLGMRES
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  KSPType, PARAMETER :: PETSC_KSPDGMRES = KSPDGMRES
  KSPType, PARAMETER :: PETSC_KSPPGMRES = KSPPGMRES
#endif
  KSPType, PARAMETER :: PETSC_KSPTCQMR = KSPTCQMR
  KSPType, PARAMETER :: PETSC_KSPBCGS = KSPBCGS
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  KSPType, PARAMETER :: PETSC_KSPIBCGS = KSPIBCGS
  KSPType, PARAMETER :: PETSC_KSPFBCGS = KSPFBCGS
  KSPType, PARAMETER :: PETSC_KSPFBCGSR = KSPFBCGSR
#endif
  KSPType, PARAMETER :: PETSC_KSPBCGSL = KSPBCGSL
  KSPType, PARAMETER :: PETSC_KSPCGS = KSPCGS
  KSPType, PARAMETER :: PETSC_KSPTFQMR = KSPTFQMR
  KSPType, PARAMETER :: PETSC_KSPCR = KSPCR
  KSPType, PARAMETER :: PETSC_KSPLSQR = KSPLSQR
  KSPType, PARAMETER :: PETSC_KSPPREONLY = KSPPREONLY
  KSPType, PARAMETER :: PETSC_KSPQCG = KSPQCG
  KSPType, PARAMETER :: PETSC_KSPBICG = KSPBICG
  KSPType, PARAMETER :: PETSC_KSPMINRES = KSPMINRES
  KSPType, PARAMETER :: PETSC_KSPSYMMLQ = KSPSYMMLQ
  KSPType, PARAMETER :: PETSC_KSPLCD = KSPLCD
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  KSPType, PARAMETER :: PETSC_KSPPYTHON = KSPPYTHON
  KSPType, PARAMETER :: PETSC_KSPGCR = KSPGCR
  KSPType, PARAMETER :: PETSC_KSPSPECEST = KSPSPECEST
#endif
  
  !MatAssembly types
  MatAssemblyType, PARAMETER :: PETSC_MAT_FLUSH_ASSEMBLY = MAT_FLUSH_ASSEMBLY
  MatAssemblyType, PARAMETER :: PETSC_MAT_FINAL_ASSEMBLY = MAT_FINAL_ASSEMBLY

  !MatDuplicate types
  MatDuplicateOption, PARAMETER :: PETSC_MAT_DO_NOT_COPY_VALUES = MAT_DO_NOT_COPY_VALUES
  MatDuplicateOption, PARAMETER :: PETSC_MAT_COPY_VALUES = MAT_COPY_VALUES
  MatDuplicateOption, PARAMETER :: PETSC_MAT_SHARE_NONZERO_PATTERN = MAT_SHARE_NONZERO_PATTERN

  !MatInfo types
  MatInfo, PARAMETER :: PETSC_MAT_INFO_SIZE = MAT_INFO_SIZE
  MatInfo, PARAMETER :: PETSC_MAT_INFO_BLOCK_SIZE = MAT_INFO_BLOCK_SIZE
  MatInfo, PARAMETER :: PETSC_MAT_INFO_NZ_ALLOCATED = MAT_INFO_NZ_ALLOCATED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_NZ_USED = MAT_INFO_NZ_USED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_NZ_UNNEEDED = MAT_INFO_NZ_UNNEEDED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_MEMORY = MAT_INFO_MEMORY
  MatInfo, PARAMETER :: PETSC_MAT_INFO_ASSEMBLIES = MAT_INFO_ASSEMBLIES
  MatInfo, PARAMETER :: PETSC_MAT_INFO_MALLOCS = MAT_INFO_MALLOCS
  MatInfo, PARAMETER :: PETSC_MAT_INFO_FILL_RATIO_GIVEN = MAT_INFO_FILL_RATIO_GIVEN
  MatInfo, PARAMETER :: PETSC_MAT_INFO_FILL_RATIO_NEEDED = MAT_INFO_FILL_RATIO_NEEDED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_FACTOR_MALLOCS = MAT_INFO_FACTOR_MALLOCS  

  !MatInfoType types
  MatInfoType, PARAMETER :: PETSC_MAT_LOCAL = MAT_LOCAL
  MatInfoType, PARAMETER :: PETSC_MAT_GLOBAL_MAX = MAT_GLOBAL_MAX
  MatInfoType, PARAMETER :: PETSC_MAT_GLOBAL_SUM = MAT_GLOBAL_SUM
  
  !MatOption types
  MatOption, PARAMETER :: PETSC_MAT_ROW_ORIENTED = MAT_ROW_ORIENTED
  MatOption, PARAMETER :: PETSC_MAT_NEW_NONZERO_LOCATIONS = MAT_NEW_NONZERO_LOCATIONS
  MatOption, PARAMETER :: PETSC_MAT_SYMMETRIC = MAT_SYMMETRIC
  MatOption, PARAMETER :: PETSC_MAT_STRUCTURALLY_SYMMETRIC = MAT_STRUCTURALLY_SYMMETRIC
  MatOption, PARAMETER :: PETSC_MAT_NEW_DIAGONALS = MAT_NEW_DIAGONALS
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_OFF_PROC_ENTRIES = MAT_IGNORE_OFF_PROC_ENTRIES
  MatOption, PARAMETER :: PETSC_MAT_NEW_NONZERO_LOCATION_ERR = MAT_NEW_NONZERO_LOCATION_ERR
  MatOption, PARAMETER :: PETSC_MAT_NEW_NONZERO_ALLOCATION_ERR = MAT_NEW_NONZERO_ALLOCATION_ERR
  MatOption, PARAMETER :: PETSC_MAT_USE_HASH_TABLE = MAT_USE_HASH_TABLE
  MatOption, PARAMETER :: PETSC_MAT_KEEP_NONZERO_PATTERN = MAT_KEEP_NONZERO_PATTERN
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_ZERO_ENTRIES = MAT_IGNORE_ZERO_ENTRIES
  MatOption, PARAMETER :: PETSC_MAT_USE_INODES = MAT_USE_INODES
  MatOption, PARAMETER :: PETSC_MAT_HERMITIAN = MAT_HERMITIAN
  MatOption, PARAMETER :: PETSC_MAT_SYMMETRY_ETERNAL = MAT_SYMMETRY_ETERNAL
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  MatOption, PARAMETER :: PETSC_MAT_DUMMY = MAT_DUMMY
#endif  
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_LOWER_TRIANGULAR = MAT_IGNORE_LOWER_TRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_ERROR_LOWER_TRIANGULAR = MAT_ERROR_LOWER_TRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_GETROW_UPPERTRIANGULAR = MAT_GETROW_UPPERTRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR = MAT_UNUSED_NONZERO_LOCATION_ERR
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR <= 4 )
  MatOption, PARAMETER :: PETSC_MAT_CHECK_COMPRESSED_ROW = MAT_CHECK_COMPRESSED_ROW
#endif  
#else
  MatOption, PARAMETER :: PETSC_MAT_USE_COMPRESSEDROW = MAT_USE_COMPRESSEDROW
#endif  
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  MatOption, PARAMETER :: PETSC_MAT_SPD = MAT_SPD
  MatOption, PARAMETER :: PETSC_MAT_NO_OFF_PROC_ENTRIES = MAT_NO_OFF_PROC_ENTRIES
  MatOption, PARAMETER :: PETSC_MAT_NO_OFF_PROC_ZERO_ROWS = MAT_NO_OFF_PROC_ZERO_ROWS
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  MatOption, PARAMETER :: PETSC_NUM_MAT_OPTIONS = MAT_OPTION_MAX
#else
  MatOption, PARAMETER :: PETSC_NUM_MAT_OPTIONS = NUM_MAT_OPTIONS
#endif
  
  !MatStructure types
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5 )
  MatStructure, PARAMETER :: PETSC_SAME_PRECONDITIONER = SAME_PRECONDITIONER
#endif  
  MatStructure, PARAMETER :: PETSC_SAME_NONZERO_PATTERN = SAME_NONZERO_PATTERN
  MatStructure, PARAMETER :: PETSC_DIFFERENT_NONZERO_PATTERN = DIFFERENT_NONZERO_PATTERN
  MatStructure, PARAMETER :: PETSC_SUBSET_NONZERO_PATTERN = SUBSET_NONZERO_PATTERN

  !MatColoring types
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_NATURAL = MATCOLORINGNATURAL
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_SL = MATCOLORINGSL
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_LF = MATCOLORINGLF
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_ID = MATCOLORINGID
#else  
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_NATURAL = MATCOLORING_NATURAL
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_SL = MATCOLORING_SL
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_LF = MATCOLORING_LF
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_ID = MATCOLORING_ID
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_GREEDY = MATCOLORINGGREEDY
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_JP = MATCOLORINGJP
#endif
  
  !Mat types
#if ( PETSC_VERSION_MAJOR == 2 )
  MatType, PARAMETER :: PETSC_AIJMUMPS = MATAIJMUMPS
#endif
  
  !PC types
  PCType, PARAMETER ::  PETSC_PCNONE = PCNONE
  PCType, PARAMETER ::  PETSC_PCJACOBI = PCJACOBI
  PCType, PARAMETER ::  PETSC_PCSOR = PCSOR
  PCType, PARAMETER ::  PETSC_PCLU = PCLU
  PCType, PARAMETER ::  PETSC_PCSHELL = PCSHELL
  PCType, PARAMETER ::  PETSC_PCBJACOBI = PCBJACOBI
  PCType, PARAMETER ::  PETSC_PCMG = PCMG
  PCType, PARAMETER ::  PETSC_PCEISENSTAT = PCEISENSTAT
  PCType, PARAMETER ::  PETSC_PCILU = PCILU
  PCType, PARAMETER ::  PETSC_PCICC = PCICC
  PCType, PARAMETER ::  PETSC_PCASM = PCASM
  PCType, PARAMETER ::  PETSC_PCKSP = PCKSP
  PCType, PARAMETER ::  PETSC_PCCOMPOSITE = PCCOMPOSITE
  PCType, PARAMETER ::  PETSC_PCREDUNDANT = PCREDUNDANT
  PCType, PARAMETER ::  PETSC_PCSPAI = PCSPAI
  PCType, PARAMETER ::  PETSC_PCNN = PCNN
  PCType, PARAMETER ::  PETSC_PCCHOLESKY = PCCHOLESKY
  PCType, PARAMETER ::  PETSC_PCPBJACOBI = PCPBJACOBI
  PCType, PARAMETER ::  PETSC_PCMAT = PCMAT
  PCType, PARAMETER ::  PETSC_PCHYPRE = PCHYPRE
  PCType, PARAMETER ::  PETSC_PCFIELDSPLIT = PCFIELDSPLIT
  PCType, PARAMETER ::  PETSC_PCML = PCML
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  PCType, PARAMETER ::  PETSC_PCGASM = PCGASM
  PCType, PARAMETER ::  PETSC_PCPARMS = PCPARMS
  PCType, PARAMETER ::  PETSC_PCTFS = PCTFS
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR <= 3 )
  PCType, PARAMETER ::  PETSC_PCPROMETHEUS = PCPROMETHEUS
#endif
  PCType, PARAMETER ::  PETSC_PCGALERKIN = PCGALERKIN
  PCType, PARAMETER ::  PETSC_PCEXOTIC = PCEXOTIC
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 4 )
  PCType, PARAMETER ::  PETSC_PCHMPI = PCHMPI
  PCType, PARAMETER ::  PETSC_PCASA = PCASA
#endif  
  PCType, PARAMETER ::  PETSC_PCSUPPORTGRAPH = PCSUPPORTGRAPH
  PCType, PARAMETER ::  PETSC_PCCP = PCCP
  PCType, PARAMETER ::  PETSC_PCVFVT = PCBFBT
  PCType, PARAMETER ::  PETSC_PCLSC = PCLSC
  PCType, PARAMETER ::  PETSC_PCPYTHON = PCPYTHON
  PCType, PARAMETER ::  PETSC_PCPFMG = PCPFMG
  PCType, PARAMETER ::  PETSC_PCSYSPFMG = PCSYSPFMG
  PCType, PARAMETER ::  PETSC_PCREDISTRIBUTE = PCREDISTRIBUTE
  PCType, PARAMETER ::  PETSC_PCSACUSP = PCSACUSP
  PCType, PARAMETER ::  PETSC_PCSACUSPPOLY = PCSACUSPPOLY
  PCType, PARAMETER ::  PETSC_PCBICGSTABCUSP = PCBICGSTABCUSP
  PCType, PARAMETER ::  PETSC_PCSVD = PCSVD
  PCType, PARAMETER ::  PETSC_PCAINVCUSP = PCAINVCUSP
#else
  PCType, PARAMETER ::  PETSC_PCMILU = PCMILU
#endif

  !Matrix Solver Package types
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR <= 3 )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SPOOLES = MATSOLVERSPOOLES
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PLAPACK = MATSOLVERPLAPACK
#endif
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU = MATSOLVERSUPERLU
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU_DIST = MATSOLVERSUPERLU_DIST
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_UMFPACK = MATSOLVERUMFPACK
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_CHOLMOD = MATSOLVERCHOLMOD
#endif 
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_ESSL = MATSOLVERESSL
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_LUSOL = MATSOLVERLUSOL
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_MUMPS = MATSOLVERMUMPS
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PASTIX = MATSOLVERPASTIX
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_MATLAB = MATSOLVERMATLAB
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PETSC = MATSOLVERPETSC
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_BAS = MATSOLVERBAS
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_CUSPARSE = MATSOLVERCUSPARSE
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_BSTRM = MATSOLVERBSTRM
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SBSTRM = MATSOLVERSBSTRM
#endif 
#else
#if ( PETSC_VERSION_MAJOR == 3 )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SPOOLES = MAT_SOLVER_SPOOLES
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU = MAT_SOLVER_SUPERLU
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU_DIST = MAT_SOLVER_SUPERLU_DIST
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_UMFPACK = MAT_SOLVER_UMFPACK
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_ESSL = MAT_SOLVER_ESSL
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_LUSOL = MAT_SOLVER_LUSOL
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_MUMPS = MAT_SOLVER_MUMPS
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_DSCPACK = MAT_SOLVER_DSCPACK
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_MATLAB = MAT_SOLVER_MATLAB
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PETSC = MAT_SOLVER_PETSC
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PLAPACK = MAT_SOLVER_PLAPACK
#if ( PETSC_VERSION_MINOR >= 1 )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PASTIX = MAT_SOLVER_PASTIX
#endif
#endif
#endif
  
  !SNES converged types
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_FNORM_ABS = SNES_CONVERGED_FNORM_ABS
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_FNORM_RELATIVE = SNES_CONVERGED_FNORM_RELATIVE
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_PNORM_RELATIVE = SNES_CONVERGED_PNORM_RELATIVE
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_SNORM_RELATIVE = SNES_CONVERGED_SNORM_RELATIVE
#endif
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_ITS = SNES_CONVERGED_ITS
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_TR_DELTA = SNES_CONVERGED_TR_DELTA
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FUNCTION_DOMAIN = SNES_DIVERGED_FUNCTION_DOMAIN
#endif
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FUNCTION_COUNT = SNES_DIVERGED_FUNCTION_COUNT
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LINEAR_SOLVE = SNES_DIVERGED_LINEAR_SOLVE
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FNORM_NAN = SNES_DIVERGED_FNORM_NAN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_MAX_IT = SNES_DIVERGED_MAX_IT
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LS_FAILURE = SNES_DIVERGED_LINE_SEARCH
#else
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LS_FAILURE = SNES_DIVERGED_LS_FAILURE
#endif
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LOCAL_MIN = SNES_DIVERGED_LOCAL_MIN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_ITERATING = SNES_CONVERGED_ITERATING
  
  !SNES types
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 4 )
  SNESType, PARAMETER :: PETSC_SNESLS = SNESNEWTONLS
  SNESType, PARAMETER :: PETSC_SNESTR = SNESNEWTONTR
#else
  SNESType, PARAMETER :: PETSC_SNESLS = SNESLS
  SNESType, PARAMETER :: PETSC_SNESTR = SNESTR
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  SNESType, PARAMETER :: PETSC_SNESPYTHON = SNESPYTHON
#endif
  SNESType, PARAMETER :: PETSC_SNESTEST = SNESTEST
  SNESType, PARAMETER :: PETSC_SNESNRICHARDSON = SNESNRICHARDSON
  SNESType, PARAMETER :: PETSC_SNESKSPONLY = SNESKSPONLY
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 4 )
  SNESType, PARAMETER :: PETSC_SNESVIRS = SNESVINEWTONRSLS
  SNESType, PARAMETER :: PETSC_SNESVISS = SNESVINEWTONSSLS
#else
  SNESType, PARAMETER :: PETSC_SNESVIRS = SNESVIRS
  SNESType, PARAMETER :: PETSC_SNESVISS = SNESVISS
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  SNESType, PARAMETER :: PETSC_SNESNGMRES = SNESNGMRES
  SNESType, PARAMETER :: PETSC_SNESQN = SNESQN
  SNESType, PARAMETER :: PETSC_SNESSHELL = SNESSHELL
  SNESType, PARAMETER :: PETSC_SNESNCG = SNESNCG
  SNESType, PARAMETER :: PETSC_SNESFAS = SNESFAS
  SNESType, PARAMETER :: PETSC_SNESMS = SNESMS
#endif

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  !SNES norm schedules
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_DEFAULT = SNES_NORM_DEFAULT
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_NONE = SNES_NORM_NONE
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_ALWAYS = SNES_NORM_ALWAYS
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_INITIAL_ONLY = SNES_NORM_INITIAL_ONLY
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_FINAL_ONLY = SNES_NORM_FINAL_ONLY
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_INITIAL_FINAL_ONLY = SNES_NORM_INITIAL_FINAL_ONLY
#else
  !SNES norm types
  SNESNormType, PARAMETER :: PETSC_SNES_NORM_DEFAULT = SNES_NORM_DEFAULT
  SNESNormType, PARAMETER :: PETSC_SNES_NORM_NONE = SNES_NORM_NONE
  SNESNormType, PARAMETER :: PETSC_SNES_NORM_FUNCTION = SNES_NORM_FUNCTION
  SNESNormType, PARAMETER :: PETSC_SNES_NORM_INITIAL_ONLY = SNES_NORM_INITIAL_ONLY
  SNESNormType, PARAMETER :: PETSC_SNES_NORM_FINAL_ONLY = SNES_NORM_FINAL_ONLY
  SNESNormType, PARAMETER :: PETSC_SNES_NORM_INITIAL_FINAL_ONLY = SNES_NORM_INITIAL_FINAL_ONLY
#endif
  
  !SNES line search type
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_BASIC = SNESLINESEARCHBASIC
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_BT = SNESLINESEARCHBT
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_L2 = SNESLINESEARCHL2
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_CP = SNESLINESEARCHCP
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_SHELL = SNESLINESEARCHSHELL
  
  !SNES line search order  
  SNESLineSearchOrder, PARAMETER :: PETSC_SNES_LINESEARCH_LINEAR = SNES_LINESEARCH_ORDER_LINEAR
  SNESLineSearchOrder, PARAMETER :: PETSC_SNES_LINESEARCH_QUADRATIC = SNES_LINESEARCH_ORDER_QUADRATIC
  SNESLineSearchOrder, PARAMETER :: PETSC_SNES_LINESEARCH_CUBIC = SNES_LINESEARCH_ORDER_CUBIC
#else
  !SNES line search types
  INTEGER(INTG), PARAMETER :: PETSC_SNES_LINESEARCH_NONORMS = 1
  INTEGER(INTG), PARAMETER :: PETSC_SNES_LINESEARCH_NO = 2
  INTEGER(INTG), PARAMETER :: PETSC_SNES_LINESEARCH_QUADRATIC = 3
  INTEGER(INTG), PARAMETER :: PETSC_SNES_LINESEARCH_CUBIC = 4  
#endif

   !TS types
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  TSType, PARAMETER :: PETSC_TS_EULER = TSEULER
  TSType, PARAMETER :: PETSC_TS_BEULER = TSBEULER
  TSType, PARAMETER :: PETSC_TS_PSEUDO = TSPSEUDO
  TSType, PARAMETER :: PETSC_TS_CRANK_NICHOLSON = TSCN
  TSType, PARAMETER :: PETSC_TS_SUNDIALS = TSSUNDIALS
  TSType, PARAMETER :: PETSC_TS_RUNGE_KUTTA = TSRK
  TSType, PARAMETER :: PETSC_TS_PYTHON = TSPYTHON
  TSType, PARAMETER :: PETSC_TS_THETA = TSTHETA
  TSType, PARAMETER :: PETSC_TS_ALPHA = TSGL
  TSType, PARAMETER :: PETSC_TS_SSP = TSSSP
  TSType, PARAMETER :: PETSC_TS_ARKIMEX = TSARKIMEX
#else
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 1 )
  TSType, PARAMETER :: PETSC_TS_EULER = TS_EULER
  TSType, PARAMETER :: PETSC_TS_BEULER = TS_BEULER
  TSType, PARAMETER :: PETSC_TS_PSEUDO = TS_PSEUDO
  TSType, PARAMETER :: PETSC_TS_SUNDIALS = TS_SUNDIALS
  TSType, PARAMETER :: PETSC_TS_CRANK_NICHOLSON = TS_CRANK_NICHOLSON
  TSType, PARAMETER :: PETSC_TS_RUNGE_KUTTA = TS_RUNGE_KUTTA
#else
  TSType, PARAMETER :: PETSC_TS_EULER = TSEULER
  TSType, PARAMETER :: PETSC_TS_BEULER = TSBEULER
  TSType, PARAMETER :: PETSC_TS_PSEUDO = TSPSEUDO
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  TSType, PARAMETER :: PETSC_TS_CRANK_NICHOLSON = TSCN
#else
  TSType, PARAMETER :: PETSC_TS_CRANK_NICHOLSON = TSCRANK_NICHOLSON
#endif  
  TSType, PARAMETER :: PETSC_TS_SUNDIALS = TSSUNDIALS
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  TSType, PARAMETER :: PETSC_TS_RUNGE_KUTTA = TSRK
#else
  TSType, PARAMETER :: PETSC_TS_RUNGE_KUTTA = TSRUNGE_KUTTA
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  TSType, PARAMETER :: PETSC_TS_PYTHON = TSPYTHON
#endif  
  TSType, PARAMETER :: PETSC_TS_THETA = TSTHETA
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  TSType, PARAMETER :: PETSC_TS_ALPHA = TSALPHA
#endif  
  TSType, PARAMETER :: PETSC_TS_GL = TSGL
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  TSType, PARAMETER :: PETSC_TS_SSP = TSSSP
  TSType, PARAMETER :: PETSC_TS_ARKIMEX = TSARKIMEX
  TSType, PARAMETER :: PETSC_TS_ROSW = TSROSW
  TSType, PARAMETER :: PETSC_TS_EIMEX = TSEIMEX
#endif  
#endif
#endif

  !TS convergence flags
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  TSConvergedReason, PARAMETER :: PETSC_TS_CONVERGED_ITERATING = TS_CONVERGED_ITERATING
  TSConvergedReason, PARAMETER :: PETSC_TS_CONVERGED_TIME = TS_CONVERGED_TIME
  TSConvergedReason, PARAMETER :: PETSC_TS_CONVERGED_ITS = TS_CONVERGED_ITS
  TSConvergedReason, PARAMETER :: PETSC_TS_DIVERGED_NONLINEAR_SOLVE = TS_DIVERGED_NONLINEAR_SOLVE
  TSConvergedReason, PARAMETER :: PETSC_TS_DIVERGED_STEP_REJECTED = TS_DIVERGED_STEP_REJECTED
#endif
  
  !TS problem types
  TSProblemType, PARAMETER :: PETSC_TS_LINEAR = TS_LINEAR
  TSProblemType, PARAMETER :: PETSC_TS_NONLINEAR = TS_NONLINEAR
  
  !TS Sundials types
  TSSundialsType, PARAMETER :: PETSC_SUNDIALS_ADAMS = SUNDIALS_ADAMS
  TSSundialsType, PARAMETER :: PETSC_SUNDIALS_BDF = SUNDIALS_BDF

  !TS Sundials Gram Schmidt Type
#if ( PETSC_VERSION_MAJOR == 3 )
  TSSundialsGramSchmidtType, PARAMETER :: PETSC_SUNDIALS_MODIFIED_GS = SUNDIALS_MODIFIED_GS
  TSSundialsGramSchmidtType, PARAMETER :: PETSC_SUNDIALS_CLASSICAL_GS = SUNDIALS_CLASSICAL_GS
#else
  TSSundialsGramSchmitdType, PARAMETER :: PETSC_SUNDIALS_MODIFIED_GS = SUNDIALS_MODIFIED_GS
  TSSundialsGramSchmitdType, PARAMETER :: PETSC_SUNDIALS_CLASSICAL_GS = SUNDIALS_CLASSICAL_GS
#endif
  
  !Module types

  !Module variables

  LOGICAL, SAVE :: PETSC_HANDLE_ERROR

  !Interfaces

  INTERFACE

    SUBROUTINE ISDestroy(indexset,ierr)
      IS indexset
      PetscInt ierr
    END SUBROUTINE ISDestroy
    
    SUBROUTINE ISColoringDestroy(iscoloring,ierr)
      ISColoring iscoloring
      PetscInt ierr
    END SUBROUTINE ISColoringDestroy
    
    SUBROUTINE ISLocalToGlobalMappingApply(ctx,type,nin,idxin,nout,idxout,ierr)
      ISLocalToGlobalMapping ctx
      ISGlobalToLocalMappingType type
      PetscInt nin
      PetscInt idxin(*)
      PetscInt nout
      PetscInt idxout(*)
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingApply
    
    SUBROUTINE ISLocalToGlobalMappingApplyIS(ctx,isin,isout,ierr)
      ISLocalToGlobalMapping ctx
      IS isin
      IS isout
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingApplyIS
    
    SUBROUTINE ISLocalToGlobalMappingCreate(comm,N,globalnum,ctx,ierr)
      MPI_Comm comm
      PetscInt N
      PetscInt globalnum(*)
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingCreate
    
    SUBROUTINE ISLocalToGlobalMappingDestroy(ctx,ierr)
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE ISLocalToGlobalMappingDestroy
    
    SUBROUTINE KSPCreate(comm,ksp,ierr)
      MPI_Comm comm
          KSP ksp
      PetscInt ierr
    END SUBROUTINE KSPCreate
    
    SUBROUTINE KSPDestroy(ksp,ierr)
      KSP ksp
      PetscInt ierr
    END SUBROUTINE KSPDestroy
    
    SUBROUTINE KSPGetConvergedReason(ksp,reason,ierr)
      KSP ksp
      KSPConvergedReason reason
      PetscInt ierr
    END SUBROUTINE KSPGetConvergedReason
    
    SUBROUTINE KSPGetIterationNumber(ksp,its,ierr)
      KSP ksp
      PetscInt its
      PetscInt ierr
    END SUBROUTINE KSPGetIterationNumber
    
    SUBROUTINE KSPGetPC(ksp,pc,ierr)
      KSP ksp
      PC pc
      PetscInt ierr
    END SUBROUTINE KSPGetPC
    
    SUBROUTINE KSPGetResidualNorm(ksp,rnorm,ierr)
      KSP ksp
      PetscReal rnorm
      PetscInt ierr
    END SUBROUTINE KSPGetResidualNorm
    
    SUBROUTINE KSPGMRESSetRestart(ksp,restart,ierr)
      KSP ksp
      PetscInt restart
      PetscInt ierr
    END SUBROUTINE KSPGMRESSetRestart
    
    SUBROUTINE KSPSetFromOptions(ksp,ierr)
      KSP ksp
      PetscInt ierr
    END SUBROUTINE KSPSetFromOptions
    
    SUBROUTINE KSPSetInitialGuessNonzero(ksp,flg,ierr)
      KSP ksp
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
      PetscBool flg
#else
      PetscTruth flg
#endif
      PetscInt ierr
    END SUBROUTINE KSPSetInitialGuessNonzero
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    SUBROUTINE KSPSetOperators(ksp,amat,pmat,ierr)
      KSP ksp
      Mat amat
      Mat pmat
      PetscInt ierr
    END SUBROUTINE KSPSetOperators
#else    
    SUBROUTINE KSPSetOperators(ksp,Amat,Pmat,flag,ierr)
      KSP ksp
      Mat Amat
      Mat Pmat
      MatStructure flag
      PetscInt ierr
    END SUBROUTINE KSPSetOperators
#endif    
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    SUBROUTINE KSPSetReusePreconditioner(ksp,flag,ierr)
      KSP ksp
      PetscBool flag
      PetscInt ierr
    END SUBROUTINE KSPSetReusePreconditioner
#endif    
    
    SUBROUTINE KSPSetTolerances(ksp,rtol,atol,dtol,maxits,ierr)
      KSP ksp
      PetscReal rtol
      PetscReal atol
      PetscReal dtol
      PetscInt maxits
      PetscInt ierr
    END SUBROUTINE KSPSetTolerances
    
    SUBROUTINE KSPSetType(ksp,method,ierr)
      KSP ksp
      KSPType method
      PetscInt ierr
    END SUBROUTINE KSPSetType
    
    SUBROUTINE KSPSetUp(ksp,ierr)
      KSP ksp
      PetscInt ierr
    END SUBROUTINE KSPSetUp
    
    SUBROUTINE KSPSolve(ksp,b,x,ierr)
      KSP ksp
      Vec b
      Vec x
      PetscInt ierr
    END SUBROUTINE KSPSolve
    
    SUBROUTINE MatAssemblyBegin(A,assemblytype,ierr)
      Mat A
      MatAssemblyType assemblytype
      PetscInt ierr
    END SUBROUTINE MatAssemblyBegin
    
    SUBROUTINE MatAssemblyEnd(A,assemblytype,ierr)
      Mat A
      MatAssemblyType assemblytype
      PetscInt ierr
    END SUBROUTINE MatAssemblyEnd
    
    SUBROUTINE MatCreate(comm,A,ierr)
      MPI_Comm comm
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreate

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
    SUBROUTINE MatCreateAIJ(comm,localm,localn,globalm,globaln,diagnumbernzperrow,diagnumbernzeachrow,offdiagnumbernzperrow, &
      & offdiagnumbernzeachrow,A,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt localn
      PetscInt globalm
      PetscInt globaln
      PetscInt diagnumbernzperrow
      PetscInt diagnumbernzeachrow(*)
      PetscInt offdiagnumbernzperrow
      PetscInt offdiagnumbernzeachrow(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateAIJ
#else
    SUBROUTINE MatCreateMPIAIJ(comm,localm,localn,globalm,globaln,diagnumbernzperrow,diagnumbernzeachrow,offdiagnumbernzperrow, &
      & offdiagnumbernzeachrow,A,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt localn
      PetscInt globalm
      PetscInt globaln
      PetscInt diagnumbernzperrow
      PetscInt diagnumbernzeachrow(*)
      PetscInt offdiagnumbernzperrow
      PetscInt offdiagnumbernzeachrow(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateMPIAIJ
#endif    
        
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
    SUBROUTINE MatCreateDense(comm,localm,localn,globalm,globaln,matrixdata,A,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt localn
      PetscInt globalm
      PetscInt globaln
      PetscScalar matrixdata(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateDense
#else
    SUBROUTINE MatCreateMPIDense(comm,localm,localn,globalm,globaln,matrixdata,A,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt localn
      PetscInt globalm
      PetscInt globaln
      PetscScalar matrixdata(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateMPIDense
#endif    
        
    SUBROUTINE MatCreateSeqAIJ(comm,m,n,numbernzperrow,numbernzeachrow,A,ierr)
      MPI_Comm comm
      PetscInt m
      PetscInt n
      PetscInt numbernzperrow
      PetscInt numbernzeachrow(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateSeqAIJ
          
    SUBROUTINE MatCreateSeqDense(comm,m,n,matrixdata,A,ierr)
      MPI_Comm comm
      PetscInt m
      PetscInt n
      PetscScalar matrixdata(*)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreateSeqDense

    SUBROUTINE MatSetType(A,matrixtype,ierr)
      Mat A
      MatType matrixtype
      PetscInt ierr
    END SUBROUTINE MatSetType

    SUBROUTINE MatDestroy(A,ierr)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatDestroy
    
    SUBROUTINE MatFDColoringCreate(A,iscoloring,fdcoloring,ierr)
      Mat A
      ISColoring iscoloring
      MatFDColoring fdcoloring
      PetscInt ierr
    END SUBROUTINE MatFDColoringCreate
    
    SUBROUTINE MatFDColoringDestroy(fdcoloring,ierr)
      MatFDColoring fdcoloring
      PetscInt ierr
    END SUBROUTINE MatFDColoringDestroy
    
    SUBROUTINE MatFDColoringSetFromOptions(fdcoloring,ierr)
      MatFDColoring fdcoloring
      PetscInt ierr
    END SUBROUTINE MatFDColoringSetFromOptions

    SUBROUTINE MatFDColoringSetParameters(fdcoloring,rerror,umin,ierr)
      MatFDColoring fdcoloring
      PetscScalar rerror
      PetscScalar umin
      PetscInt ierr
    END SUBROUTINE MatFDColoringSetParameters
    
    SUBROUTINE MatFDColoringSetFunction(fdcoloring,ffunction,ctx,ierr)
      USE TYPES
      MatFDColoring fdcoloring
      EXTERNAL ffunction
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE MatFDColoringSetFunction
    
    SUBROUTINE MatFDColoringSetFunctionSNES(fdcoloring,ffunction,ctx,ierr)
      USE TYPES
      MatFDColoring fdcoloring
      EXTERNAL ffunction
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE MatFDColoringSetFunctionSNES
    
    SUBROUTINE MatGetArray(A,mat_data,mat_offset,ierr)
      Mat A
      PetscScalar mat_data(1)
      PetscOffset mat_offset
      PetscInt ierr
    END SUBROUTINE MatGetArray
    
    SUBROUTINE MatGetArrayF90(A,mat_data,ierr)
      Mat A
      PetscScalar, POINTER :: mat_data(:,:)
      PetscInt ierr
    END SUBROUTINE MatGetArrayF90
    
    SUBROUTINE MatGetColoring(A,coloring_type,iscoloring,ierr)
      Mat A
      MatColoringType coloring_type
      ISColoring iscoloring 
      PetscInt ierr
    END SUBROUTINE MatGetColoring
    
    SUBROUTINE MatGetInfo(A,flag,info,ierr)
      Mat A
      MatInfoType flag
      MatInfo info(*)
      PetscInt ierr
    END SUBROUTINE MatGetInfo
    
    SUBROUTINE MatGetOwnershipRange(A,firstrow,lastrow,ierr)
      Mat A
      PetscInt firstrow
      PetscInt lastrow
      PetscInt ierr
    END SUBROUTINE MatGetOwnershipRange
    
    SUBROUTINE MatGetRow(A,row,ncols,cols,values,ierr)
      Mat A
      PetscInt row
      PetscInt ncols
      PetscInt cols(*)
      PetscScalar values(*)
      PetscInt ierr
    END SUBROUTINE MatGetRow
    
    SUBROUTINE MatGetValues(A,m,idxm,n,idxn,values,ierr)
      Mat A
      PetscInt m
      PetscInt idxm(*)
      PetscInt n
      PetscInt idxn(*)
      PetscScalar values(*)
      PetscInt ierr
    END SUBROUTINE MatGetValues
    
    SUBROUTINE MatRestoreArray(A,mat_data,mat_offset,ierr)
      Mat A
      PetscScalar mat_data(1)
      PetscOffset mat_offset
      PetscInt ierr
    END SUBROUTINE MatRestoreArray
        
    SUBROUTINE MatRestoreArrayF90(A,mat_data,ierr)
      Mat A
      PetscScalar, POINTER :: mat_data(:,:)
      PetscInt ierr
    END SUBROUTINE MatRestoreArrayF90
    
   SUBROUTINE MatRestoreRow(A,row,ncols,cols,values,ierr)
      Mat A
      PetscInt row
      PetscInt ncols
      PetscInt cols(*)
      PetscScalar values(*)
      PetscInt ierr
    END SUBROUTINE MatRestoreRow
    
    SUBROUTINE MatSetLocalToGlobalMapping(A,ctx,ierr)
      Mat A
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE MatSetLocalToGlobalMapping

    SUBROUTINE MatSetOption(A,option,flag,ierr)
      Mat A
      MatOption option
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
      PetscBool flag
#else
      PetscTruth flag
#endif
      PetscInt ierr
    END SUBROUTINE MatSetOption

    SUBROUTINE MatSetSizes(A,localm,localn,globalM,globalN,ierr)
      Mat A
      PetscInt localm
      PetscInt localn
      PetscInt globalM
      PetscInt globalN
      PetscInt ierr
    END SUBROUTINE MatSetSizes
    
    SUBROUTINE MatSetValue(A,row,col,value,insertmode,ierr)
      Mat A
      PetscInt row
      PetscInt col
      PetscScalar value
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE MatSetValue
    
    SUBROUTINE MatSetValues(A,m,mindices,n,nindices,values,insertmode,ierr)
      Mat A
      PetscInt m
      PetscInt mindices(*)
      PetscInt n
      PetscInt nindices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE MatSetValues
    
    SUBROUTINE MatSetValuesLocal(A,m,mindices,n,nindices,values,insertmode,ierr)
      Mat A
      PetscInt m
      PetscInt mindices(*)
      PetscInt n
      PetscInt nindices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE MatSetValuesLocal
    
    SUBROUTINE MatSetValueLocal(A,row,col,value,insertmode,ierr)
      Mat A
      PetscInt row
      PetscInt col
      PetscScalar value
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE MatSetValueLocal
    
    SUBROUTINE MatView(A,v,ierr)
      Mat A
      PetscViewer v
      PetscInt ierr
    END SUBROUTINE MatView

    SUBROUTINE MatZeroEntries(A,ierr)
      Mat A
      PetscInt ierr
    END SUBROUTINE MatZeroEntries

#if ( PETSC_MAJOR_VERSION >= 3 ) 
    SUBROUTINE PCFactorSetMatSolverPackage(pc,solverpackage,ierr)
      PC pc
      MatSolverPackage solverpackage
      PetscInt ierr
    END SUBROUTINE PCFactorSetMatSolverPackage

    SUBROUTINE PCFactorSetUpMatSolverPackage(pc,ierr)
      PC pc
      PetscInt ierr
    END SUBROUTINE PCFactorSetUpMatSolverPackage

    SUBROUTINE PCFactorGetMatrix(pc,A,ierr)
      PC pc
      Mat A
      PetscInt ierr
    END SUBROUTINE PCFactorGetMatrix

    SUBROUTINE MatMumpsSetIcntl(A,icntl,ival,ierr)
      Mat A
      PetscInt icntl
      PetscInt ival
      PetscInt ierr
    END SUBROUTINE MatMumpsSetIcntl

#if ( PETSC_VERSION_MINOR >= 4 )
    SUBROUTINE MatMumpsSetCntl(A,icntl,val,ierr)
      Mat A
      PetscInt icntl
      PetscReal val
      PetscInt ierr
    END SUBROUTINE MatMumpsSetCntl
#endif
#endif
    
    SUBROUTINE PCSetType(pc,method,ierr)
      PC pc
      PCType method
      PetscInt ierr
    END SUBROUTINE PCSetType
    
    SUBROUTINE PetscFinalize(ierr)
      PetscInt ierr
    END SUBROUTINE PetscFinalize
    
    SUBROUTINE PetscInitialize(file,ierr)
      PetscChar(*) file
      PetscInt ierr
    END SUBROUTINE PetscInitialize

    SUBROUTINE PetscPopSignalHandler(ierr)
      PetscInt ierr
    END SUBROUTINE PetscPopSignalHandler

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
    SUBROUTINE PetscLogView(viewer,ierr)
      PetscViewer viewer
      PetscInt ierr
    END SUBROUTINE PetscLogView
#else
    SUBROUTINE PetscLogPrintSummary(comm,file,ierr)
      MPI_Comm comm
      PetscChar(*) file
      PetscInt ierr
    END SUBROUTINE PetscLogPrintSummary
#endif

    SUBROUTINE SNESCreate(comm,snes,ierr)
      MPI_Comm comm
      SNES snes
      PetscInt ierr
    END SUBROUTINE SNESCreate

    SUBROUTINE SNESDestroy(snes,ierr)
      SNES snes
      PetscInt ierr
    END SUBROUTINE SNESDestroy

    SUBROUTINE SNESSetConvergenceTest(snes,cfunction,ctx,destroy,ierr)
      USE TYPES
      SNES snes
      EXTERNAL cfunction
      TYPE(SOLVER_TYPE), POINTER :: ctx
      EXTERNAL destroy
      PetscInt ierr
    END SUBROUTINE SNESSetConvergenceTest

    SUBROUTINE SNESGetConvergedReason(snes,reason,ierr)
      SNES snes
      SNESConvergedReason reason
      PetscInt ierr
    END SUBROUTINE SNESGetConvergedReason
    
    SUBROUTINE SNESGetFunctionNorm(snes,fnorm,ierr)
      SNES snes
      PetscReal fnorm
      PetscInt ierr
    END SUBROUTINE SNESGetFunctionNorm

    SUBROUTINE SNESGetSolutionUpdate(snes,solutionUpdate,ierr)
      SNES snes
      Vec solutionUpdate
      PetscInt ierr
    END SUBROUTINE SNESGetSolutionUpdate

    SUBROUTINE SNESSetFunctionNorm(snes,fnorm,ierr)
      SNES snes
      PetscReal fnorm
      PetscInt ierr
    END SUBROUTINE SNESSetFunctionNorm

    SUBROUTINE SnesLineSearchSetNorms(snes,xnorm,fnorm,ynorm,ierr)
      SNES snes
      PetscReal xnorm
      PetscReal fnorm
      PetscReal ynorm
      PetscInt ierr
    END SUBROUTINE SnesLineSearchSetNorms

    SUBROUTINE SnesLineSearchGetNorms(linesearch,xnorm,fnorm,ynorm,ierr)
      SNESLineSearch linesearch
      PetscReal xnorm
      PetscReal fnorm
      PetscReal ynorm
      PetscInt ierr
    END SUBROUTINE SnesLineSearchGetNorms

    SUBROUTINE SNESGetIterationNumber(snes,iter,ierr)
      SNES snes
      PetscInt iter
      PetscInt ierr
    END SUBROUTINE SNESGetIterationNumber

    SUBROUTINE SNESGetKSP(snes,ksp,ierr)
      SNES snes
      KSP ksp
      PetscInt ierr
    END SUBROUTINE SNESGetKSP

    SUBROUTINE SNESLineSearchSet(snes,func,lsctx,ierr)
      SNES snes
      EXTERNAL func
      PetscFortranAddr lsctx
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSet
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
    SUBROUTINE SnesLineSearchSetMonitor(linesearch,flag,ierr)
      SNESLineSearch linesearch
      PetscBool flag
      PetscInt ierr
    END SUBROUTINE SnesLineSearchSetMonitor
#endif

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
    SUBROUTINE SNESLineSearchSetComputeNorms(linesearch,flag,ierr)
      SNESLineSearch linesearch
      PetscBool flag
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSetComputeNorms

    SUBROUTINE SnesLineSearchComputeNorms(linesearch,ierr)
      SNESLineSearch linesearch
      PetscInt ierr
    END SUBROUTINE SnesLineSearchComputeNorms

    SUBROUTINE SNESLineSearchSetOrder(linesearch,linesearchorder,ierr)
      SNESLineSearch linesearch
      SNESLineSearchOrder linesearchorder
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSetOrder

    SUBROUTINE SNESLineSearchBTSetAlpha(linesearch,alpha,ierr)
      SNESLineSearch linesearch
      PetscReal alpha
      PetscInt ierr
    END SUBROUTINE SNESLineSearchBTSetAlpha

    SUBROUTINE SNESLineSearchSetTolerances(linesearch,steptol,maxstep,rtol,atol,ltol,maxIt,ierr)
      SNESLineSearch linesearch
      PetscReal steptol
      PetscReal maxstep
      PetscReal rtol
      PetscReal atol
      PetscReal ltol
      PetscInt maxIt
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSetTolerances
#endif
    
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
    SUBROUTINE SNESLineSearchSetParams(snes,alpha,maxstep,steptol,ierr)
      SNES snes
      PetscReal alpha
      PetscReal maxstep
      PetscReal steptol
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSetParams
#endif    
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
#if ( PETSC_VERSION_MINOR >= 4 )
    SUBROUTINE SNESGetLineSearch(snes,linesearch,ierr)
#else
    SUBROUTINE SNESGetSNESLineSearch(snes,linesearch,ierr)
#endif
      SNES snes
      SNESLineSearch linesearch
      PetscInt ierr
#if ( PETSC_VERSION_MINOR >= 4 )
    END SUBROUTINE SNESGetLineSearch
#else
    END SUBROUTINE SNESGetSNESLineSearch
#endif

    SUBROUTINE SNESLineSearchSetType(linesearch,linesearchtype,ierr)
      SNESLineSearch linesearch
      SNESLineSearchType linesearchtype
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSetType
#endif
    
    SUBROUTINE SNESMonitorSet(snes,mfunction,mctx,monitordestroy,ierr)
      USE TYPES
      SNES snes
      EXTERNAL mfunction
      TYPE(SOLVER_TYPE), POINTER :: mctx
      EXTERNAL monitordestroy
      PetscInt ierr
    END SUBROUTINE SNESMonitorSet

     SUBROUTINE SNESSetFromOptions(snes,ierr)
      SNES snes
      PetscInt ierr
    END SUBROUTINE SNESSetFromOptions

    SUBROUTINE SNESSetFunction(snes,f,ffunction,ctx,ierr)
      USE TYPES
      SNES snes
      Vec f
      EXTERNAL ffunction
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE SNESSetFunction

    SUBROUTINE SNESGetJacobian(snes,A,B,Jfunction,ctx,ierr)
      USE TYPES
      SNES snes
      Mat A
      Mat B      
      EXTERNAL Jfunction
!       ISLocalToGlobalMapping ctx
!       TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ctx
      PetscInt ierr
    END SUBROUTINE SNESGetJacobian

    SUBROUTINE SNESSetJacobian(snes,A,B,Jfunction,ctx,ierr)
      USE TYPES
      SNES snes
      Mat A
      Mat B      
      EXTERNAL Jfunction
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE SNESSetJacobian

    SUBROUTINE SNESSetKSP(snes,ksp,ierr)
      SNES snes
      KSP ksp
      PetscInt ierr
    END SUBROUTINE SNESSetKSP

    SUBROUTINE SNESSetTolerances(snes,abstol,rtol,stol,maxit,maxf,ierr)
      SNES snes
      PetscReal abstol
      PetscReal rtol
      PetscReal stol
      PetscInt maxit
      PetscInt maxf
      PetscInt ierr
    END SUBROUTINE SNESSetTolerances
    
    SUBROUTINE SNESSetTrustRegionTolerance(snes,trtol,ierr)
      SNES snes
      PetscReal trtol
      PetscInt ierr
    END SUBROUTINE SNESSetTrustRegionTolerance

    SUBROUTINE SNESSetType(snes,method,ierr)
      SNES snes
      SNESType method
      PetscInt ierr
    END SUBROUTINE SNESSetType

    SUBROUTINE SNESLineSearchGetVecs(linesearch,x,f,y,w,g,ierr)
      SNESLineSearch linesearch
      Vec x
      Vec f
      Vec y
      Vec w
      Vec g
      PetscInt ierr
    END SUBROUTINE SNESLineSearchGetVecs

#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 4 )
    SUBROUTINE SNESSetNormType(snes,normtype,ierr)
      SNES snes
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
      SNESNormSchedule normtype
#else
      SNESNormType normtype
#endif
      PetscInt ierr
    END SUBROUTINE SNESSetNormType
#endif

    SUBROUTINE SNESSolve(snes,b,x,ierr)
      SNES snes
      Vec b
      Vec x
      PetscInt ierr
    END SUBROUTINE SNESSolve

    SUBROUTINE TSCreate(comm,ts,ierr)
      MPI_Comm comm
      TS ts
      PetscInt ierr
    END SUBROUTINE TSCreate

    SUBROUTINE TSDestroy(ts,ierr)
      TS ts
      PetscInt ierr
    END SUBROUTINE TSDestroy

    SUBROUTINE TSGetApplicationContext(ts,userP,ierr)
      USE TYPES
      TS ts
      TYPE(SOLVER_TYPE), POINTER :: userP
      PetscInt ierr
    END SUBROUTINE TSGetApplicationContext

    SUBROUTINE TSMonitorSet(ts,mfunction,mctx,monitordestroy,ierr)
      USE TYPES
      TS ts
      EXTERNAL mfunction
      TYPE(SOLVER_TYPE), POINTER :: mctx
      EXTERNAL monitordestroy
      PetscInt ierr
    END SUBROUTINE TSMonitorSet

    SUBROUTINE TSSetApplicationContext(ts,userP,ierr)
      USE TYPES
      TS ts
      TYPE(SOLVER_TYPE), POINTER :: userP
      PetscInt ierr
    END SUBROUTINE TSSetApplicationContext

    SUBROUTINE TSSetDuration(ts,maxsteps,maxtime,ierr)
      TS ts
      PetscInt maxsteps
      PetscReal maxtime
      PetscInt ierr
    END SUBROUTINE TSSetDuration

    SUBROUTINE TSSetFromOptions(ts,ierr)
      TS ts
      PetscInt ierr
    END SUBROUTINE TSSetFromOptions

    SUBROUTINE TSSetInitialTimeStep(ts,initial_time,time_step,ierr)
      TS ts
      PetscReal initial_time
      PetscReal time_step
      PetscInt ierr
    END SUBROUTINE TSSetInitialTimeStep
    
    SUBROUTINE TSSetProblemType(ts,probtype,ierr)
      TS ts
      TSProblemType probtype
      PetscInt ierr
    END SUBROUTINE TSSetProblemType
    
    SUBROUTINE TSSetRHSFunction(ts,rhsfunc,ctx,ierr)
      USE TYPES
      TS ts
      EXTERNAL rhsfunc
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE TSSetRHSFunction
    
    SUBROUTINE TSSetTimeStep(ts,time_step,ierr)
      TS ts
      PetscReal time_step
      PetscInt ierr
    END SUBROUTINE TSSetTimeStep
    
    SUBROUTINE TSSetType(ts,tstype,ierr)
      TS ts
      TSType tstype
      PetscInt ierr
    END SUBROUTINE TSSetType
    
    SUBROUTINE TSSolve(ts,x,ierr)
      TS ts
      Vec x
      PetscInt ierr
    END SUBROUTINE TSSolve

    SUBROUTINE TSStep(ts,steps,ptime,ierr)
      TS ts
      PetscInt steps
      PetscReal ptime
      PetscInt ierr
    END SUBROUTINE TSStep

    SUBROUTINE VecAssemblyBegin(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecAssemblyBegin

    SUBROUTINE VecAssemblyEnd(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecAssemblyEnd

    SUBROUTINE VecCopy(x,y,ierr)
      Vec x
      Vec y
      PetscInt ierr
    END SUBROUTINE VecCopy

    SUBROUTINE VecCreate(comm,x,ierr)
      MPI_Comm comm
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreate

    SUBROUTINE VecCreateGhost(comm,localm,globalm,nghost,ghosts,x,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt globalm
      PetscInt nghost
      PetscInt ghosts(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateGhost

    SUBROUTINE VecCreateGhostWithArray(comm,localm,globalm,nghost,ghosts,array,x,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt globalm
      PetscInt nghost
      PetscInt ghosts(*)
      PetscScalar array(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateGhostWithArray

    SUBROUTINE VecCreateMPI(comm,localm,globalm,x,ierr)
      MPI_Comm comm
      PetscInt localm
      PetscInt globalm
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateMPI

    SUBROUTINE VecCreateMPIWithArray(comm,localn,globaln,array,x,ierr)
      MPI_Comm comm
      PetscInt localn
      PetscInt globaln
      PetscScalar array(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateMPIWithArray

    SUBROUTINE VecCreateSeq(comm,m,x,ierr)
      MPI_Comm comm
      PetscInt m
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateSeq

    SUBROUTINE VecCreateSeqWithArray(comm,n,array,x,ierr)
      MPI_Comm comm
      PetscInt n
      PetscScalar array(*)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecCreateSeqWithArray

    SUBROUTINE VecDestroy(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecDestroy

    SUBROUTINE VecDuplicate(old,new,ierr)
      Vec old,new
      PetscInt ierr
    END SUBROUTINE VecDuplicate

    SUBROUTINE VecDot(x,y,val,ierr)
      Vec x
      Vec y
      PetscScalar val
    END SUBROUTINE VecDot

    SUBROUTINE VecGetArray(x,vec_data,vec_offset,ierr)
      Vec x
      PetscScalar vec_data(1)
      PetscOffset vec_offset
      PetscInt ierr
    END SUBROUTINE VecGetArray
    
    SUBROUTINE VecGetArrayF90(x,vec_data,ierr)
      Vec x
      PetscScalar, POINTER :: vec_data(:)
      PetscInt ierr
    END SUBROUTINE VecGetArrayF90

    SUBROUTINE VecGetLocalSize(x,size,ierr)
      Vec x
      PetscInt size
      PetscInt ierr
    END SUBROUTINE VecGetLocalSize

    SUBROUTINE VecGetOwnershipRange(x,low,high,ierr)
      Vec x
      PetscInt low
      PetscInt high
      PetscInt ierr
    END SUBROUTINE VecGetOwnershipRange

    SUBROUTINE VecGetSize(x,size,ierr)
      Vec x
      PetscInt size
      PetscInt ierr
    END SUBROUTINE VecGetSize

    SUBROUTINE VecGetValues(x,n,indices,values,ierr)
      Vec x
      PetscInt n
      PetscInt indices(*)
      PetscScalar values(*)
      PetscInt ierr
    END SUBROUTINE VecGetValues

    SUBROUTINE VecGhostGetLocalForm(g,l,ierr)
      Vec g
      Vec l
      PetscInt ierr
    END SUBROUTINE VecGhostGetLocalForm

    SUBROUTINE VecGhostRestoreLocalForm(g,l,ierr)
      Vec g
      Vec l
      PetscInt ierr
    END SUBROUTINE VecGhostRestoreLocalForm

   SUBROUTINE VecGhostUpdateBegin(x,insertmode,scattermode,ierr)
      Vec x
      InsertMode insertmode
      ScatterMode scattermode
      PetscInt ierr
    END SUBROUTINE VecGhostUpdateBegin

    SUBROUTINE VecGhostUpdateEnd(x,insertmode,scattermode,ierr)
      Vec x
      InsertMode insertmode
      ScatterMode scattermode
      PetscInt ierr
    END SUBROUTINE VecGhostUpdateEnd

    SUBROUTINE VecRestoreArray(x,vec_data,vec_offset,ierr)
      Vec x
      PetscScalar vec_data(1)
      PetscOffset vec_offset
      PetscInt ierr
    END SUBROUTINE VecRestoreArray

    SUBROUTINE VecRestoreArrayF90(x,vec_data,ierr)
      Vec x
      PetscScalar, POINTER :: vec_data(:)
      PetscInt ierr
    END SUBROUTINE VecRestoreArrayF90

    SUBROUTINE VecScale(x,alpha,ierr)
      Vec x
      PetscScalar alpha
      PetscInt ierr
    END SUBROUTINE VecScale

    SUBROUTINE VecSet(x,value,ierr)
      Vec x
      PetscScalar value
      PetscInt ierr
    END SUBROUTINE VecSet

    SUBROUTINE VecSetFromOptions(x,ierr)
      Vec x
      PetscInt ierr
    END SUBROUTINE VecSetFromOptions

    SUBROUTINE VecSetLocalToGlobalMapping(v,ctx,ierr)
      Vec v
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE VecSetLocalToGlobalMapping

    SUBROUTINE VecSetSizes(x,localm,globalm,ierr)
      Vec x
      PetscInt localm,globalm
      PetscInt ierr
    END SUBROUTINE VecSetSizes

    SUBROUTINE VecSetValues(x,n,indices,values,insertmode,ierr)
      Vec x
      PetscInt n
      PetscInt indices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE VecSetValues

    SUBROUTINE VecSetValuesLocal(x,n,indices,values,insertmode,ierr)
      Vec x
      PetscInt n
      PetscInt indices(*)
      PetscScalar values(*)
      InsertMode insertmode
      PetscInt ierr
    END SUBROUTINE VecSetValuesLocal

    SUBROUTINE VecView(x,v,ierr)
      Vec x
      PetscViewer v
      PetscInt ierr
    END SUBROUTINE VecView

  END INTERFACE

  INTERFACE PETSC_SNESSETJACOBIAN
    MODULE PROCEDURE PETSC_SNESSETJACOBIAN_SOLVER  
    MODULE PROCEDURE PETSC_SNESSETJACOBIAN_MATFDCOLORING
  END INTERFACE !PETSC_SNESSETJACOBIAN

  INTERFACE PETSC_SNESGETJACOBIAN
    MODULE PROCEDURE PETSC_SNESGETJACOBIAN_SOLVER  
    MODULE PROCEDURE PETSC_SNESGETJACOBIAN_SPECIAL
  END INTERFACE !PETSC_SNESSETJACOBIAN

  PUBLIC PETSC_TRUE,PETSC_FALSE
  
  PUBLIC PETSC_NULL_CHARACTER,PETSC_NULL_INTEGER,PETSC_NULL_DOUBLE,PETSC_NULL_OBJECT, &
    & PETSC_NULL_FUNCTION,PETSC_NULL_SCALAR,PETSC_NULL_REAL
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 4 )
  PUBLIC PETSC_NULL
#endif
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
  PUBLIC PETSC_NULL_TRUTH
#else
  PUBLIC PETSC_NULL_BOOL
#endif

  PUBLIC PETSC_ADD_VALUES,PETSC_INSERT_VALUES,PETSC_COMM_WORLD,PETSC_COMM_SELF,PETSC_DECIDE, &
    & PETSC_SCATTER_FORWARD,PETSC_SCATTER_REVERSE

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  PUBLIC PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_REAL
#else
  PUBLIC PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_DOUBLE_PRECISION
#endif  
  
  PUBLIC PETSC_KSPRICHARDSON,PETSC_KSPCG,PETSC_KSPCGNE,PETSC_KSPSTCG,PETSC_KSPGMRES,PETSC_KSPFGMRES, &
    & PETSC_KSPLGMRES,PETSC_KSPTCQMR,PETSC_KSPBCGS,PETSC_KSPBCGSL,PETSC_KSPCGS,PETSC_KSPTFQMR,PETSC_KSPCR,PETSC_KSPLSQR, &
    & PETSC_KSPPREONLY,PETSC_KSPQCG,PETSC_KSPBICG,PETSC_KSPMINRES,PETSC_KSPSYMMLQ,PETSC_KSPLCD

#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
  PUBLIC PETSC_KSPCHEBYCHEV
#else
  PUBLIC PETSC_KSPCHEBYSHEV
#endif

#if ( PETSC_VERSION_MAJOR == 2 )
  PUBLIC PETSC_AIJMUMPS
#endif
  
  PUBLIC PETSC_PCNONE,PETSC_PCJACOBI,PETSC_PCSOR,PETSC_PCLU,PETSC_PCSHELL,PETSC_PCBJACOBI,PETSC_PCMG,PETSC_PCEISENSTAT, &
    & PETSC_PCILU,PETSC_PCICC,PETSC_PCASM,PETSC_PCKSP,PETSC_PCCOMPOSITE,PETSC_PCREDUNDANT,PETSC_PCSPAI, &
    & PETSC_PCNN,PETSC_PCCHOLESKY,PETSC_PCPBJACOBI,PETSC_PCMAT,PETSC_PCHYPRE,PETSC_PCFIELDSPLIT,PETSC_PCML
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR < 4 )
  PUBLIC PETSC_PCPROMETHEUS
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR <= 4 )
  PUBLIC PETSC_PCHMPI,PETSC_PCASA
#endif  
  PUBLIC PETSC_PCGASM,PETSC_PCPARMS,PETSC_PCTFS,PETSC_PCGALERKIN,PETSC_PCEXOTIC, &
    & PETSC_PCSUPPORTGRAPH,PETSC_PCCP,PETSC_PCVFVT,PETSC_PCLSC,PETSC_PCPYTHON,PETSC_PCPFMG,PETSC_PCSYSPFMG, &
    & PETSC_PCREDISTRIBUTE,PETSC_PCSACUSP,PETSC_PCSACUSPPOLY,PETSC_PCBICGSTABCUSP,PETSC_PCSVD,PETSC_PCAINVCUSP
#else
  PUBLIC PETSC_PCMILU
#endif
  
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 4 )
  PUBLIC PETSC_SAME_PRECONDITIONER
#endif
  
  PUBLIC PETSC_SAME_NONZERO_PATTERN,PETSC_DIFFERENT_NONZERO_PATTERN,PETSC_SUBSET_NONZERO_PATTERN

  PUBLIC PETSC_ISINITIALISE,PETSC_ISFINALISE,PETSC_ISDESTROY

  PUBLIC PETSC_ISCOLORINGINITIALISE,PETSC_ISCOLORINGFINALISE,PETSC_ISCOLORINGDESTROY
  
  PUBLIC PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE,PETSC_ISLOCALTOGLOBALMAPPINGFINALISE,PETSC_ISLOCALTOGLOBALMAPPINGAPPLY, &
    & PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS,PETSC_ISLOCALTOGLOBALMAPPINGCREATE,PETSC_ISLOCALTOGLOBALMAPPINGDESTROY

  PUBLIC PETSC_KSP_CONVERGED_RTOL,PETSC_KSP_CONVERGED_ATOL,PETSC_KSP_CONVERGED_ITS,PETSC_KSP_CONVERGED_CG_NEG_CURVE, &
    & PETSC_KSP_CONVERGED_CG_CONSTRAINED,PETSC_KSP_CONVERGED_STEP_LENGTH,PETSC_KSP_CONVERGED_HAPPY_BREAKDOWN, &
    & PETSC_KSP_DIVERGED_NULL,PETSC_KSP_DIVERGED_ITS,PETSC_KSP_DIVERGED_DTOL,PETSC_KSP_DIVERGED_BREAKDOWN, &
    & PETSC_KSP_DIVERGED_BREAKDOWN_BICG,PETSC_KSP_DIVERGED_NONSYMMETRIC,PETSC_KSP_DIVERGED_INDEFINITE_PC, &
    & PETSC_KSP_DIVERGED_NAN,PETSC_KSP_DIVERGED_INDEFINITE_MAT,PETSC_KSP_CONVERGED_ITERATING
  
  PUBLIC PETSC_KSPCREATE,PETSC_KSPDESTROY,PETSC_KSPGETCONVERGEDREASON,PETSC_KSPGETITERATIONNUMBER,PETSC_KSPGETPC, &
    & PETSC_KSPGETRESIDUALNORM,PETSC_KSPGMRESSETRESTART,PETSC_KSPFINALISE,PETSC_KSPINITIALISE,PETSC_KSPSETFROMOPTIONS, &
    & PETSC_KSPSETINITIALGUESSNONZERO,PETSC_KSPSETOPERATORS,PETSC_KSPSETTYPE,PETSC_KSPSETUP,PETSC_KSPSETTOLERANCES, &
    & PETSC_KSPSOLVE
  
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  PUBLIC PETSc_KspSetReusePreconditioner
#endif
  
#if ( PETSC_VERSION_MAJOR < 3 )
  PUBLIC MAT_COLUMN_ORIENTED,MAT_COLUMNS_SORTED,MAT_ROWS_SORTED,MAT_FINAL_ASSEMBLY,MAT_FLUSH_ASSEMBLY, &
    & MAT_NO_NEW_NONZERO_LOCATIONS
#endif

  PUBLIC PETSC_MAT_FLUSH_ASSEMBLY,PETSC_MAT_FINAL_ASSEMBLY

  PUBLIC PETSC_MAT_DO_NOT_COPY_VALUES,PETSC_MAT_COPY_VALUES,PETSC_MAT_SHARE_NONZERO_PATTERN
  
  PUBLIC PETSC_MAT_INFO_SIZE,PETSC_MAT_INFO_BLOCK_SIZE,PETSC_MAT_INFO_NZ_ALLOCATED,PETSC_MAT_INFO_NZ_USED, &
    & PETSC_MAT_INFO_NZ_UNNEEDED,PETSC_MAT_INFO_MEMORY,PETSC_MAT_INFO_ASSEMBLIES,PETSC_MAT_INFO_MALLOCS, &
    & PETSC_MAT_INFO_FILL_RATIO_GIVEN,PETSC_MAT_INFO_FILL_RATIO_NEEDED,PETSC_MAT_INFO_FACTOR_MALLOCS

  PUBLIC PETSC_MAT_LOCAL,PETSC_MAT_GLOBAL_MAX,PETSC_MAT_GLOBAL_SUM

  PUBLIC PETSC_MAT_ROW_ORIENTED,PETSC_MAT_NEW_NONZERO_LOCATIONS,PETSC_MAT_SYMMETRIC,PETSC_MAT_STRUCTURALLY_SYMMETRIC, &
    & PETSC_MAT_NEW_DIAGONALS,PETSC_MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_MAT_NEW_NONZERO_LOCATION_ERR, &
    & PETSC_MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_MAT_USE_HASH_TABLE,PETSC_MAT_KEEP_NONZERO_PATTERN, &
    & PETSC_MAT_IGNORE_ZERO_ENTRIES,PETSC_MAT_USE_INODES,PETSC_MAT_HERMITIAN,PETSC_MAT_SYMMETRY_ETERNAL, &
    & PETSC_MAT_IGNORE_LOWER_TRIANGULAR,PETSC_MAT_GETROW_UPPERTRIANGULAR, &
    & PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR,PETSC_NUM_MAT_OPTIONS
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR <= 4 )  
  PUBLIC PETSC_MAT_CHECK_COMPRESSED_ROW
#endif  
  PUBLIC PETSC_MAT_SPD,PETSC_MAT_NO_OFF_PROC_ENTRIES,PETSC_MAT_NO_OFF_PROC_ZERO_ROWS
#else
  PUBLIC PETSC_MAT_USE_COMPRESSEDROW
#endif
  
  PUBLIC PETSC_MATCOLORING_NATURAL,PETSC_MATCOLORING_SL,PETSC_MATCOLORING_LF,PETSC_MATCOLORING_ID

  PUBLIC PETSC_MATINITIALISE,PETSC_MATFINALISE,PETSC_MATASSEMBLYBEGIN,PETSC_MATASSEMBLYEND,PETSC_MATCREATE, &
    & PETSC_MATCREATESEQAIJ,PETSC_MATCREATESEQDENSE,PETSC_MATDESTROY, &
    & PETSC_MATFDCOLORINGCREATE,PETSC_MATFDCOLORINGDESTROY,PETSC_MATFDCOLORINGFINALISE,PETSC_MATFDCOLORINGINITIALISE, &
    & PETSC_MATFDCOLORINGSETFROMOPTIONS,PETSC_MATFDCOLORINGSETFUNCTION,PETSC_MATFDCOLORINGSETFUNCTIONSNES, &
    & PETSC_MATFDCOLORINGSETPARAMETERS, &
    & PETSC_MATGETARRAY,PETSC_MATGETARRAYF90,PETSC_MATGETINFO,PETSC_MATGETCOLORING,PETSC_MATGETOWNERSHIPRANGE, &
    & PETSC_MATGETROW,PETSC_MATGETVALUES,PETSC_MATRESTOREARRAY,PETSC_MATRESTOREARRAYF90,PETSC_MATRESTOREROW, &
    & PETSC_MATSETLOCALTOGLOBALMAPPING,PETSC_MATSETOPTION,PETSC_MATSETSIZES,PETSC_MATSETVALUE,PETSC_MATSETVALUES, &
    & PETSC_MATSETVALUELOCAL,PETSC_MATSETVALUESLOCAL,PETSC_MATVIEW,PETSC_MATZEROENTRIES,PETSC_MATSETTYPE
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  PUBLIC PETSC_MATCREATEAIJ,PETSC_MATCREATEDENSE
#else
  PUBLIC PETSC_MATCREATEMPIAIJ,PETSC_MATCREATEMPIDENSE
#endif
  
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR < 4 )
  PUBLIC PETSC_MAT_SOLVER_SPOOLES,PETSC_MAT_SOLVER_PLAPACK
#endif
  PUBLIC PETSC_MAT_SOLVER_SUPERLU,PETSC_MAT_SOLVER_SUPERLU_DIST,PETSC_MAT_SOLVER_UMFPACK, &
    & PETSC_MAT_SOLVER_ESSL,PETSC_MAT_SOLVER_LUSOL,PETSC_MAT_SOLVER_MUMPS,PETSC_MAT_SOLVER_MATLAB,PETSC_MAT_SOLVER_PASTIX, &
    & PETSC_MAT_SOLVER_PETSC
#else
#if ( PETSC_VERSION_MAJOR == 3 )
  PUBLIC PETSC_MAT_SOLVER_SPOOLES,PETSC_MAT_SOLVER_SUPERLU,PETSC_MAT_SOLVER_SUPERLU_DIST,PETSC_MAT_SOLVER_UMFPACK, &
    & PETSC_MAT_SOLVER_ESSL,PETSC_MAT_SOLVER_LUSOL,PETSC_MAT_SOLVER_MUMPS,PETSC_MAT_SOLVER_DSCPACK,PETSC_MAT_SOLVER_MATLAB, &
    & PETSC_MAT_SOLVER_PETSC,PETSC_MAT_SOLVER_PLAPACK
#if ( PETSC_VERSION_MINOR >= 1 )
  PUBLIC PETSC_MAT_SOLVER_PASTIX
#endif
#endif
#endif
  
  PUBLIC PETSC_PCINITIALISE,PETSC_PCFINALISE,PETSC_PCSETTYPE

#if ( PETSC_VERSION_MAJOR == 3 )
  PUBLIC PETSC_PCFACTORSETMATSOLVERPACKAGE
  PUBLIC Petsc_PCFactorSetUpMatSolverPackage
  PUBLIC Petsc_PCFactorGetMatrix
  PUBLIC Petsc_MatMumpsSetIcntl
#if ( PETSC_VERSION_MINOR >= 4 )
  PUBLIC Petsc_MatMumpsSetCntl
#endif
#endif
  
  PUBLIC PETSC_TS_EULER,PETSC_TS_BEULER,PETSC_TS_PSEUDO,PETSC_TS_SUNDIALS,PETSC_TS_CRANK_NICHOLSON,PETSC_TS_RUNGE_KUTTA
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  PUBLIC PETSC_TS_PYTHON,PETSC_TS_THETA,PETSC_TS_ALPHA,PETSC_TS_SSP,PETSC_TS_ARKIMEX
#else
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 1 )
  PUBLIC PETSC_TS_THETA,PETSC_TS_GL
#endif
#endif

  PUBLIC PETSC_TS_LINEAR,PETSC_TS_NONLINEAR

  PUBLIC PETSC_SUNDIALS_ADAMS,PETSC_SUNDIALS_BDF,PETSC_SUNDIALS_MODIFIED_GS,PETSC_SUNDIALS_CLASSICAL_GS

  PUBLIC PETSC_TSCREATE,PETSC_TSDESTROY,PETSC_TSFINALISE,PETSC_TSINITIALISE,PETSC_TSMONITORSET, &
    & PETSC_TSSETDURATION,PETSC_TSSETFROMOPTIONS,PETSC_TSSETINITIALTIMESTEP, &
    & PETSC_TSSETPROBLEMTYPE,PETSC_TSSETRHSFUNCTION,PETSC_TSSETTIMESTEP,PETSC_TSSETTYPE,PETSC_TSSOLVE,PETSC_TSSTEP
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2 )
  PUBLIC PETSC_TSSETMATRICES
#endif
  
  PUBLIC PETSC_ERRORHANDLING_SET_OFF,PETSC_ERRORHANDLING_SET_ON
  
  PUBLIC PETSC_FINALIZE,PETSC_INITIALIZE
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  PUBLIC PETSC_LOGVIEW
#else
  PUBLIC PETSC_LOGPRINTSUMMARY
#endif
  
  PUBLIC PETSC_SNESLS,PETSC_SNESTR,PETSC_SNESTEST
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  PUBLIC PETSC_SNESPYTHON
#endif
  
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  PUBLIC PETSC_SNES_LINESEARCH_BASIC,PETSC_SNES_LINESEARCH_BT,PETSC_SNES_LINESEARCH_L2,PETSC_SNES_LINESEARCH_CP, &
    & PETSC_SNES_LINESEARCH_SHELL
  
  PUBLIC PETSC_SNES_LINESEARCH_LINEAR,PETSC_SNES_LINESEARCH_QUADRATIC,PETSC_SNES_LINESEARCH_CUBIC
#else
  PUBLIC PETSC_SNES_LINESEARCH_NONORMS,PETSC_SNES_LINESEARCH_NO,PETSC_SNES_LINESEARCH_QUADRATIC,PETSC_SNES_LINESEARCH_CUBIC
#endif

#  
  PUBLIC PETSC_SNES_NORM_DEFAULT,PETSC_SNES_NORM_NONE,PETSC_SNES_NORM_INITIAL_ONLY, &
    & PETSC_SNES_NORM_FINAL_ONLY,PETSC_SNES_NORM_INITIAL_FINAL_ONLY
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  PUBLIC PETSC_SNES_NORM_ALWAYS
#else
  PUBLIC PETSC_SNES_NORM_FUNCTION
#endif
  
  PUBLIC PETSC_SNES_CONVERGED_FNORM_ABS,PETSC_SNES_CONVERGED_FNORM_RELATIVE, &
    & PETSC_SNES_CONVERGED_ITS,PETSC_SNES_CONVERGED_TR_DELTA,PETSC_SNES_DIVERGED_FUNCTION_COUNT,PETSC_SNES_DIVERGED_LINEAR_SOLVE, &
    & PETSC_SNES_DIVERGED_FNORM_NAN,PETSC_SNES_DIVERGED_MAX_IT,PETSC_SNES_DIVERGED_LOCAL_MIN,PETSC_SNES_CONVERGED_ITERATING
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  PUBLIC PETSC_SNES_DIVERGED_FUNCTION_DOMAIN,PETSC_SNES_DIVERGED_LS_FAILURE
#else
  PUBLIC PETSC_SNES_DIVERGED_LS_FAILURE
#endif
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 2 )
  PUBLIC PETSC_SNES_CONVERGED_PNORM_RELATIVE
#endif
  
  PUBLIC PETSC_SNESFINALISE,PETSC_SNESINITIALISE,PETSC_SNESCREATE,PETSC_SNESDESTROY,PETSC_SNESGETCONVERGEDREASON, &
    & PETSC_SNESGETFUNCTIONNORM,PETSC_SNESSETFUNCTIONNORM,petsc_sneslinesearchsetnorms,petsc_sneslinesearchgetnorms, &
    & PETSC_SNESGETITERATIONNUMBER,PETSC_SNESGETKSP,PETSC_SNESMONITORSET,PETSC_SNESSETFROMOPTIONS,PETSC_SNESSETFUNCTION, &
    & PETSC_SNESSETJACOBIAN,PETSC_SNESSETTOLERANCES,PETSC_SNESSETTRUSTREGIONTOLERANCE,PETSC_SNESSETTYPE,PETSC_SNESSOLVE, &
    & PETSC_SNESSETKSP,PETSC_SNESGETJACOBIAN,PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR,PETSC_SNESDEFAULTCOMPUTEJACOBIAN, &
    & PETSC_SNESSETCONVERGENCETEST,Petsc_SnesLineSearchGetVecs,Petsc_SnesGetSolutionUpdate
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  PUBLIC Petsc_SnesLineSearchSetMonitor
#endif
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 4 )
  PUBLIC PETSC_SNESSETNORMTYPE
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  PUBLIC Petsc_SnesLineSearchFinalise,Petsc_SnesLineSearchInitialise
  PUBLIC Petsc_SnesGetSnesLineSearch,Petsc_SnesLineSearchSetComputeNorms,Petsc_SnesLineSearchComputeNorms, &
    & Petsc_SnesLineSearchSetOrder,Petsc_SnesLineSearchSetType, &
    & Petsc_SnesLineSearchBTSetAlpha,Petsc_SnesLineSearchSetTolerances
#else
  PUBLIC PETSC_SNESLINESEARCHSET,PETSC_SNESLINESEARCHSETPARAMS
#endif
  
  PUBLIC PETSC_VECINITIALISE,PETSC_VECFINALISE,PETSC_VECASSEMBLYBEGIN,PETSC_VECASSEMBLYEND,PETSC_VECCOPY,PETSC_VECCREATE, &
    & PETSC_VECCREATEGHOST,PETSC_VECCREATEGHOSTWITHARRAY,PETSC_VECCREATEMPI,PETSC_VECCREATEMPIWITHARRAY,PETSC_VECCREATESEQ, &
    & PETSC_VECCREATESEQWITHARRAY,PETSC_VECDESTROY,PETSC_VECDUPLICATE,PETSC_VECGETARRAYF90, &
    & PETSC_VECGETLOCALSIZE,PETSC_VECGETOWNERSHIPRANGE,PETSC_VECGETSIZE,PETSC_VECGETVALUES,PETSC_VECGHOSTGETLOCALFORM, &
    & PETSC_VECGHOSTRESTORELOCALFORM,PETSC_VECGHOSTUPDATEBEGIN,PETSC_VECGHOSTUPDATEEND, &
    & PETSC_VECRESTOREARRAYF90,PETSC_VECSCALE,PETSC_VECSET,PETSC_VECSETFROMOPTIONS,PETSC_VECSETLOCALTOGLOBALMAPPING, &
    & PETSC_VECSETSIZES,PETSC_VECSETVALUES,PETSC_VECSETVALUESLOCAL,PETSC_VECVIEW,Petsc_VecDot

  PUBLIC PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_DRAW_WORLD,PETSC_VIEWER_DRAW_SELF

CONTAINS

  !
  !================================================================================================================================
  !

  !>Set PETSc error handling on
  SUBROUTINE PETSC_ERRORHANDLING_SET_OFF(ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ERRORHANDLING_SET_OFF",ERR,ERROR,*999)

    PETSC_HANDLE_ERROR=.FALSE.
    
    CALL EXITS("PETSC_ERRORHANDLING_SET_OFF")
    RETURN
999 CALL ERRORS("PETSC_ERRORHANDLING_SET_OFF",ERR,ERROR)
    CALL EXITS("PETSC_ERRORHANDLING_SET_OFF")
    RETURN 1
  END SUBROUTINE PETSC_ERRORHANDLING_SET_OFF
    
  !
  !================================================================================================================================
  !

  !>Set PETSc error handling on
  SUBROUTINE PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ERRORHANDLING_SET_ON",ERR,ERROR,*999)

    PETSC_HANDLE_ERROR=.TRUE.
    
    CALL EXITS("PETSC_ERRORHANDLING_SET_ON")
    RETURN
999 CALL ERRORS("PETSC_ERRORHANDLING_SET_ON",ERR,ERROR)
    CALL EXITS("PETSC_ERRORHANDLING_SET_ON")
    RETURN 1
  END SUBROUTINE PETSC_ERRORHANDLING_SET_ON
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscFinalize routine
  SUBROUTINE PETSC_FINALIZE(ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_FINALIZE",ERR,ERROR,*999)

    CALL PetscFinalize(ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscFinalize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_FINALIZE")
    RETURN
999 CALL ERRORS("PETSC_FINALIZE",ERR,ERROR)
    CALL EXITS("PETSC_FINALIZE")
    RETURN 1
  END SUBROUTINE PETSC_FINALIZE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscInitialize routine.
  SUBROUTINE PETSC_INITIALIZE(FILE,ERR,ERROR,*)

    !Argument Variables
    CHARACTER(LEN=*), INTENT(IN) :: FILE !<Filename for PETSc options file
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_INITIALIZE",ERR,ERROR,*999)

    PETSC_HANDLE_ERROR=.TRUE.
    CALL PetscInitialize(FILE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscInitialize",ERR,ERROR,*999)
    ENDIF
    ! Disable PETSc's signal handler as we have our own OpenCMISS signal handlers in cmiss_c.c
    CALL PetscPopSignalHandler(ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscPopSignalHandler",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_INITIALIZE")
    RETURN
999 CALL ERRORS("PETSC_INITIALIZE",ERR,ERROR)
    CALL EXITS("PETSC_INITIALIZE")
    RETURN 1
  END SUBROUTINE PETSC_INITIALIZE
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc IS structure and destroy the IS
  SUBROUTINE PETSC_ISFINALISE(IS_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_IS_TYPE), INTENT(INOUT) :: IS_ !<The IS to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISFINALISE",ERR,ERROR,*999)

    IF(IS_%IS_/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_ISDESTROY(IS_,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISFINALISE")
    RETURN
999 CALL ERRORS("PETSC_ISFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_ISFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_ISFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc IS structure
  SUBROUTINE PETSC_ISINITIALISE(IS_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_IS_TYPE), INTENT(INOUT) :: IS_ !<The IS to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISINITIALISE",ERR,ERROR,*999)

    IS_%IS_=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_ISINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_ISINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_ISINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_ISINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISDestroy routine.
  SUBROUTINE PETSC_ISDESTROY(IS_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_IS_TYPE), INTENT(INOUT) :: IS_ !<The index set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISDESTROY",ERR,ERROR,*999)

    CALL ISDestroy(IS_%IS_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISDestroy",ERR,ERROR,*999)
    ENDIF
    IS_%IS_=PETSC_NULL_INTEGER
        
    CALL EXITS("PETSC_ISDESTROY")
    RETURN
999 CALL ERRORS("PETSC_ISDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_ISDESTROY")
    RETURN 1
    
  END SUBROUTINE PETSC_ISDESTROY
    
  !
  !
  !================================================================================================================================
  !

  !Finalise the PETSc ISColoring structure and destroy the ISColoring
  SUBROUTINE PETSC_ISCOLORINGFINALISE(ISCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISCOLORING_TYPE), INTENT(INOUT) :: ISCOLORING !<The ISColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISCOLORINGFINALISE",ERR,ERROR,*999)

    IF(ISCOLORING%ISCOLORING/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_ISCOLORINGDESTROY(ISCOLORING,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISCOLORINGFINALISE")
    RETURN
999 CALL ERRORS("PETSC_ISCOLORINGFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_ISCOLORINGFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_ISCOLORINGFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc ISColoring structure
  SUBROUTINE PETSC_ISCOLORINGINITIALISE(ISCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISCOLORING_TYPE), INTENT(INOUT) :: ISCOLORING !<The ISColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISCOLORINGINITIALISE",ERR,ERROR,*999)

    ISCOLORING%ISCOLORING=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_ISCOLORINGINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_ISCOLORINGINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_ISCOLORINGINITIALISE")
    RETURN 1
    
  END SUBROUTINE PETSC_ISCOLORINGINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISColoringDestroy routine.
  SUBROUTINE PETSC_ISCOLORINGDESTROY(ISCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISCOLORING_TYPE), INTENT(INOUT) :: ISCOLORING !<The index set coloring
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISCOLORINGDESTROY",ERR,ERROR,*999)

    CALL ISColoringDestroy(ISCOLORING%ISCOLORING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISColoringDestroy",ERR,ERROR,*999)
    ENDIF
    ISCOLORING%ISCOLORING=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_ISCOLORINGDESTROY")
    RETURN
999 CALL ERRORS("PETSC_ISCOLORINGDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_ISCOLORINGDESTROY")
    RETURN 1
    
  END SUBROUTINE PETSC_ISCOLORINGDESTROY
  
  !  
  !================================================================================================================================
  !

  !Finalise the PETSc ISLocalToGlobalMapping structure and destroy the ISLocalToGlobalMapping
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGFINALISE(ISLOCALTOGLOBALMAPPING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(INOUT) :: ISLOCALTOGLOBALMAPPING !<The ISLocalToGlobalMapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE",ERR,ERROR,*999)

    IF(ISLOCALTOGLOBALMAPPING%ISLOCALTOGLOBALMAPPING/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_ISLOCALTOGLOBALMAPPINGDESTROY(ISLOCALTOGLOBALMAPPING,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc ISLocalToGlobalMapping structure
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE(ISLOCALTOGLOBALMAPPING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(INOUT) :: ISLOCALTOGLOBALMAPPING !<The ISLocalToGlobalMapping to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE",ERR,ERROR,*999)

    ISLOCALTOGLOBALMAPPING%ISLOCALTOGLOBALMAPPING=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingApply routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLY(CTX,TYPE,NIN,IDXIN,NOUT,IDXOUT,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(INOUT) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(IN) :: TYPE !<The type of local to global mapping
    INTEGER(INTG), INTENT(IN) :: NIN !<The number of local indicies
    INTEGER(INTG), INTENT(IN) :: IDXIN(*)
    INTEGER(INTG), INTENT(OUT) :: NOUT
    INTEGER(INTG), INTENT(OUT) :: IDXOUT(*)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingApply(CTX%ISLOCALTOGLOBALMAPPING,TYPE,NIN,IDXIN,NOUT,IDXOUT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingApply",ERR,ERROR,*999)
    ENDIF
    CTX%ISLOCALTOGLOBALMAPPING=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY")
    RETURN 1
    
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingApplyIS routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS(CTX,ISIN,ISOUT,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(IN) :: CTX !<The local to global mapping context
    TYPE(PETSC_IS_TYPE), INTENT(IN) :: ISIN
    TYPE(PETSC_IS_TYPE), INTENT(OUT) :: ISOUT
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingApplyIS(CTX%ISLOCALTOGLOBALMAPPING,ISIN%IS_,ISOUT%IS_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingApplyIS",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingCreate routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGCREATE(COMMUNICATOR,N,GLOBALNUM,CTX,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: N !<The number of local indices
    INTEGER(INTG), INTENT(IN) :: GLOBALNUM(*) !<The global number for each local index
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(INOUT) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingCreate(COMMUNICATOR,N,GLOBALNUM,CTX%ISLOCALTOGLOBALMAPPING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingDestroy routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGDESTROY(CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(INOUT) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY",ERR,ERROR,*999)

    CALL ISLocalToGlobalMappingDestroy(CTX%ISLOCALTOGLOBALMAPPING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in ISLocalToGlobalMappingDestroy",ERR,ERROR,*999)
    ENDIF
    CTX%ISLOCALTOGLOBALMAPPING=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY")
    RETURN
999 CALL ERRORS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPCreate routine
  SUBROUTINE PETSC_KSPCREATE(COMMUNICATOR,KSP_,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator for the KSP creation
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<On exit, the Ksp information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPCREATE",ERR,ERROR,*999)

    CALL KSPCreate(COMMUNICATOR,KSP_%KSP_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPCREATE")
    RETURN
999 CALL ERRORS("PETSC_KSPCREATE",ERR,ERROR)
    CALL EXITS("PETSC_KSPCREATE")
    RETURN 1
  END SUBROUTINE PETSC_KSPCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPDestroy routine
  SUBROUTINE PETSC_KSPDESTROY(KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPDESTROY",ERR,ERROR,*999)

    CALL KSPDestroy(KSP_%KSP_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPDestroy",ERR,ERROR,*999)
    ENDIF
    KSP_%KSP_=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_KSPDESTROY")
    RETURN
999 CALL ERRORS("PETSC_KSPDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_KSPDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_KSPDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetConvergedReason routine
  SUBROUTINE PETSC_KSPGETCONVERGEDREASON(KSP_,REASON,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The KSP information
    INTEGER(INTG), INTENT(OUT) :: REASON !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPGETCONVERGEDREASON",ERR,ERROR,*999)

    CALL KSPGetConvergedReason(KSP_%KSP_,REASON,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPGetConvergedReason",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPGETCONVERGEDREASON")
    RETURN
999 CALL ERRORS("PETSC_KSPGETCONVERGEDREASON",ERR,ERROR)
    CALL EXITS("PETSC_KSPGETCONVERGEDREASON")
    RETURN 1
  END SUBROUTINE PETSC_KSPGETCONVERGEDREASON
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetIterationNumber routine
  SUBROUTINE PETSC_KSPGETITERATIONNUMBER(KSP_,ITERATION_NUMBER,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The KSP information
    INTEGER(INTG), INTENT(OUT) :: ITERATION_NUMBER !<On exit, the number of iterations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPGETITERATIONNUMBER",ERR,ERROR,*999)

    CALL KSPGetIterationNumber(KSP_%KSP_,ITERATION_NUMBER,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPGetIterationNumber",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPGETITERATIONNUMBER")
    RETURN
999 CALL ERRORS("PETSC_KSPGETITERATIONNUMBER",ERR,ERROR)
    CALL EXITS("PETSC_KSPGETITERATIONNUMBER")
    RETURN 1
  END SUBROUTINE PETSC_KSPGETITERATIONNUMBER
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetPC routine
  SUBROUTINE PETSC_KSPGETPC(KSP_,PC_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to get the PC for
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<On exit, the PC associated with the Ksp
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPGETPC",ERR,ERROR,*999)

    CALL KSPGetPC(KSP_%KSP_,PC_%PC_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPGetPC",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPGETPC")
    RETURN
999 CALL ERRORS("PETSC_KSPGETPC",ERR,ERROR)
    CALL EXITS("PETSC_KSPGETPC")
    RETURN 1
  END SUBROUTINE PETSC_KSPGETPC
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetResidualNorm routine
  SUBROUTINE PETSC_KSPGETRESIDUALNORM(KSP_,RESIDUAL_NORM,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to get the PC for
    REAL(DP), INTENT(OUT) :: RESIDUAL_NORM !<On exit, the residual norm for the KSP
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPGETRESIDUALNORM",ERR,ERROR,*999)

    CALL KSPGetResidualNorm(KSP_%KSP_,RESIDUAL_NORM,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPGetResidualNorm.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPGETRESIDUALNORM")
    RETURN
999 CALL ERRORS("PETSC_KSPGETRESIDUALNORM",ERR,ERROR)
    CALL EXITS("PETSC_KSPGETRESIDUALNORM")
    RETURN 1
  END SUBROUTINE PETSC_KSPGETRESIDUALNORM
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc KSP structure and destroy the KSP
  SUBROUTINE PETSC_KSPFINALISE(KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPFINALISE",ERR,ERROR,*999)

    IF(KSP_%KSP_/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_KSPDESTROY(KSP_,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPFINALISE")
    RETURN
999 CALL ERRORS("PETSC_KSPFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_KSPFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_KSPFINALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGMRESSetRestart routine
  SUBROUTINE PETSC_KSPGMRESSETRESTART(KSP_,RESTART,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the GMRES restart for
    INTEGER(INTG), INTENT(OUT) :: RESTART !<The restart value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPGMRESSETRESTART",ERR,ERROR,*999)

    CALL KSPGMRESSetRestart(KSP_%KSP_,RESTART,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPGMRESSetRestart.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPGMRESSETRESTART")
    RETURN
999 CALL ERRORS("PETSC_KSPGMRESSETRESTART",ERR,ERROR)
    CALL EXITS("PETSC_KSPGMRESSETRESTART")
    RETURN 1
  END SUBROUTINE PETSC_KSPGMRESSETRESTART
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc KSP structure
  SUBROUTINE PETSC_KSPINITIALISE(KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPINITIALISE",ERR,ERROR,*999)

    KSP_%KSP_=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_KSPINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_KSPINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_KSPINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_KSPINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetFromOptions routine
  SUBROUTINE PETSC_KSPSETFROMOPTIONS(KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the options for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSETFROMOPTIONS",ERR,ERROR,*999)

    CALL KSPSetFromOptions(KSP_%KSP_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSetFromOptions",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSETFROMOPTIONS")
    RETURN
999 CALL ERRORS("PETSC_KSPSETFROMOPTIONS",ERR,ERROR)
    CALL EXITS("PETSC_KSPSETFROMOPTIONS")
    RETURN 1
  END SUBROUTINE PETSC_KSPSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetInitialGuessNonzero routine
  SUBROUTINE PETSC_KSPSETINITIALGUESSNONZERO(KSP_,FLAG,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the initial guess non zero for
    LOGICAL, INTENT(IN) :: FLAG
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSETINITIALGUESSNONZERO",ERR,ERROR,*999)

    IF(FLAG) THEN
      CALL KSPSetInitialGuessNonzero(KSP_%KSP_,PETSC_TRUE,ERR)
    ELSE
      CALL KSPSetInitialGuessNonzero(KSP_%KSP_,PETSC_FALSE,ERR)
    ENDIF
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSetInitialGuessNonzero",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSETINITIALGUESSNONZERO")
    RETURN
999 CALL ERRORS("PETSC_KSPSETINITIALGUESSNONZERO",ERR,ERROR)
    CALL EXITS("PETSC_KSPSETINITIALGUESSNONZERO")
    RETURN 1
  END SUBROUTINE PETSC_KSPSETINITIALGUESSNONZERO
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  !>Buffer routine to the PETSc KSPSetOperators routine
  SUBROUTINE PETSC_KSPSetOperators(ksp_,amat,pmat,err,error,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: ksp_ !<The Ksp to set the operators for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: amat !<The matrix associated with the linear system
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: pmat !<The matrix to be used in constructing the preconditioner
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("PETSC_KSPSetOperators",err,error,*999)

    CALL KSPSetOperators(ksp_%ksp_,amat%mat,pmat%mat,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FlagError("PETSc error in KSPSetFromOperators.",err,error,*999)
    ENDIF
    
    CALL Exits("PETSC_KSPSetOperators")
    RETURN
999 CALL Errors("PETSC_KSPSetOperators",err,error)
    CALL Exits("PETSC_KSPSetOperators")
    RETURN 1
  END SUBROUTINE PETSC_KSPSetOperators
  
#else
  
  !>Buffer routine to the PETSc KSPSetOperators routine
  SUBROUTINE PETSC_KSPSETOPERATORS(KSP_,AMAT,PMAT,FLAG,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the operators for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: AMAT !<The matrix associated with the linear system
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: PMAT !<The matrix to be used in constructing the preconditioner
    MatStructure, INTENT(IN) :: FLAG !<Preconditioner matrix structure flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSETOPERATORS",ERR,ERROR,*999)

    CALL KSPSetOperators(KSP_%KSP_,AMAT%MAT,PMAT%MAT,FLAG,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSetFromOperators",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSETOPERATORS")
    RETURN
999 CALL ERRORS("PETSC_KSPSETOPERATORS",ERR,ERROR)
    CALL EXITS("PETSC_KSPSETOPERATORS")
    RETURN 1
  END SUBROUTINE PETSC_KSPSETOPERATORS
#endif
  
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  !>Buffer routine to the PETSc KSPSetReusePreconditioner routine
  SUBROUTINE PETSc_KSPSetReusePreconditioner(KSP_,flag,err,error,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the options for
    LOGICAL, INTENT(IN) :: flag
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("",err,error,*999)

    IF(flag) THEN
      CALL KSPSetReusePreconditioner(KSP_%KSP_,PETSC_TRUE,err)
    ELSE
      CALL KSPSetReusePreconditioner(KSP_%KSP_,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetReusePreconditioner.",err,error,*999)
    ENDIF
    
    CALL Exits("")
    RETURN
999 CALL Errors("",err,error)
    CALL Exits("")
    RETURN 1
  END SUBROUTINE PETSc_KSPSetReusePreconditioner
#endif
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetTolerances routine
  SUBROUTINE PETSC_KSPSETTOLERANCES(KSP_,RTOL,ATOL,DTOL,MAXITS,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the tolerances for
    REAL(DP), INTENT(IN) :: RTOL !<The relative tolerance to set
    REAL(DP), INTENT(IN) :: ATOL !<The absolution tolerance to set
    REAL(DP), INTENT(IN) :: DTOL !<The divergence tolerance to set
    INTEGER(INTG), INTENT(IN) :: MAXITS !<The maximum number of iterations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSETTOLERANCES",ERR,ERROR,*999)

    CALL KSPSetTolerances(KSP_%KSP_,RTOL,ATOL,DTOL,MAXITS,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSetTolerances",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSETTOLERANCES")
    RETURN
999 CALL ERRORS("PETSC_KSPSETTOLERANCES",ERR,ERROR)
    CALL EXITS("PETSC_KSPSETTOLERANCES")
    RETURN 1
  END SUBROUTINE PETSC_KSPSETTOLERANCES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetType routine
  SUBROUTINE PETSC_KSPSETTYPE(KSP_,METHOD,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set the type for
    KSPType, INTENT(IN) :: METHOD !<The Ksp method
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSETTYPE",ERR,ERROR,*999)

    CALL KSPSetType(KSP_%KSP_,METHOD,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSetType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSETTYPE")
    RETURN
999 CALL ERRORS("PETSC_KSPSETTYPE",ERR,ERROR)
    CALL EXITS("PETSC_KSPSETTYPE")
    RETURN 1
  END SUBROUTINE PETSC_KSPSETTYPE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetUp routine
  SUBROUTINE PETSC_KSPSETUP(KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set up
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSETUP",ERR,ERROR,*999)

    CALL KSPSetUp(KSP_%KSP_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSetUp",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSETUP")
    RETURN
999 CALL ERRORS("PETSC_KSPSETUP",ERR,ERROR)
    CALL EXITS("PETSC_KSPSETUP")
    RETURN 1
  END SUBROUTINE PETSC_KSPSETUP
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSolve routine
  SUBROUTINE PETSC_KSPSOLVE(KSP_,B,X,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The Ksp to set up
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: B !<The RHS vector
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT)  :: X !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_KSPSOLVE",ERR,ERROR,*999)

    CALL KSPSolve(KSP_%KSP_,B%VEC,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in KSPSolve",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_KSPSOLVE")
    RETURN
999 CALL ERRORS("PETSC_KSPSOLVE",ERR,ERROR)
    CALL EXITS("PETSC_KSPSOLVE")
    RETURN 1
  END SUBROUTINE PETSC_KSPSOLVE
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscLogView routine.
  SUBROUTINE PETSC_LOGVIEW(VIEWER,ERR,ERROR,*)

    !Argument Variables
    PetscViewer, INTENT(IN) :: VIEWER !<The viewer to print the log to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_LOGVIEW",ERR,ERROR,*999)

    CALL PetscLogView(VIEWER,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscLogView.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_LOGVIEW")
    RETURN
999 CALL ERRORS("PETSC_LOGVIEW",ERR,ERROR)
    CALL EXITS("PETSC_LOGVIEW")
    RETURN 1
  END SUBROUTINE PETSC_LOGVIEW
  
#else
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscLogPrintSummary routine.
  SUBROUTINE PETSC_LOGPRINTSUMMARY(COMMUNICATOR,FILE,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    CHARACTER(LEN=*), INTENT(IN) :: FILE !<Filename for the log summary
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_LOGPRINTSUMMARY",ERR,ERROR,*999)

    CALL PetscLogPrintSummary(COMMUNICATOR,FILE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PetscLogPrintSummary",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_LOGPRINTSUMMARY")
    RETURN
999 CALL ERRORS("PETSC_LOGPRINTSUMMARY",ERR,ERROR)
    CALL EXITS("PETSC_LOGPRINTSUMMARY")
    RETURN 1
  END SUBROUTINE PETSC_LOGPRINTSUMMARY
    
  !
  !================================================================================================================================
  !
  
#endif
  
  !>Buffer routine to the PETSc MatAssemblyBegin routine.
  SUBROUTINE PETSC_MATASSEMBLYBEGIN(A,ASSEMBLY_TYPE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !The matrix to assemble
    MatAssemblyType, INTENT(IN) :: ASSEMBLY_TYPE !<The assembly type 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATASSEMBLYBEGIN",ERR,ERROR,*999)

    CALL MatAssemblyBegin(A%MAT,ASSEMBLY_TYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatAssemblyBegin",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATASSEMBLYBEGIN")
    RETURN
999 CALL ERRORS("PETSC_MATASSEMBLYBEGIN",ERR,ERROR)
    CALL EXITS("PETSC_MATASSEMBLYBEGIN")
    RETURN 1
  END SUBROUTINE PETSC_MATASSEMBLYBEGIN
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatAssemblyEnd routine.
  SUBROUTINE PETSC_MATASSEMBLYEND(A,ASSEMBLY_TYPE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to assemble
    MatAssemblyType, INTENT(IN) :: ASSEMBLY_TYPE !<The assembly type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATASSEMBLYEND",ERR,ERROR,*999)

    CALL MatAssemblyEnd(A%MAT,ASSEMBLY_TYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatAssemblyEnd",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATASSEMBLYEND")
    RETURN
999 CALL ERRORS("PETSC_MATASSEMBLYEND",ERR,ERROR)
    CALL EXITS("PETSC_MATASSEMBLYEND")
    RETURN 1
  END SUBROUTINE PETSC_MATASSEMBLYEND
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreate routine.
  SUBROUTINE PETSC_MATCREATE(COMMUNICATOR,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI Communicator
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATE",ERR,ERROR,*999)

    CALL MatCreate(COMMUNICATOR,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATE    
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !>Buffer routine to the PETSc MatCreateAIJ routine.
  SUBROUTINE PETSC_MATCREATEAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
    & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<The number of global rows
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: DIAG_NUMBER_NZ_PERROW !<The maximum number of non-zeros per row in the diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: DIAG_NUMBER_NZ_EACHROW(*) !<The number of non-zeros per row in the diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: OFFDIAG_NUMBER_NZ_PERROW !<The maximum number of non-zeros per row in the off-diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: OFFDIAG_NUMBER_NZ_EACHROW(*) !<The number of non-zeros per row in the off-diagonal part of the matrix
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATEAIJ",ERR,ERROR,*999)

    CALL MatCreateAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
      & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateAIJ",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATEAIJ")
    RETURN
999 CALL ERRORS("PETSC_MATCREATEAIJ",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATEAIJ")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEAIJ
#else  
  !>Buffer routine to the PETSc MatCreateMPIAIJ routine.
  SUBROUTINE PETSC_MATCREATEMPIAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
    & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<The number of global rows
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: DIAG_NUMBER_NZ_PERROW !<The maximum number of non-zeros per row in the diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: DIAG_NUMBER_NZ_EACHROW(*) !<The number of non-zeros per row in the diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: OFFDIAG_NUMBER_NZ_PERROW !<The maximum number of non-zeros per row in the off-diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: OFFDIAG_NUMBER_NZ_EACHROW(*) !<The number of non-zeros per row in the off-diagonal part of the matrix
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATEMPIAIJ",ERR,ERROR,*999)

    CALL MatCreateMPIAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
      & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateMPIAIJ",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATEMPIAIJ")
    RETURN
999 CALL ERRORS("PETSC_MATCREATEMPIAIJ",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATEMPIAIJ")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEMPIAIJ
#endif  
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !>Buffer routine to the PETSc MatCreateDense routine.
  SUBROUTINE PETSC_MATCREATEDENSE(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<The number of global rows
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*) !<Optional, the allocated matrix data.
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATEDENSE",ERR,ERROR,*999)

    CALL MatCreateDense(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateDense",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATEDENSE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATEDENSE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATEDENSE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEDENSE
#else  
  !>Buffer routine to the PETSc MatCreateMPIDense routine.
  SUBROUTINE PETSC_MATCREATEMPIDENSE(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<The number of global rows
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*) !<Optional, the allocated matrix data.
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATEMPIDENSE",ERR,ERROR,*999)

    CALL MatCreateMPIDense(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateMPIDense",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATEMPIDENSE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATEMPIDENSE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATEMPIDENSE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEMPIDENSE
#endif
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateSeqAIJ routine.
  SUBROUTINE PETSC_MATCREATESEQAIJ(COMMUNICATOR,M,N,NUMBER_NZ_PERROW,NUMBER_NZ_EACHROW,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: M !<The number of rows
    INTEGER(INTG), INTENT(IN) :: N !<The number of columns
    INTEGER(INTG), INTENT(IN) :: NUMBER_NZ_PERROW !<The maximum number of non-zeros per row
    INTEGER(INTG), INTENT(IN) :: NUMBER_NZ_EACHROW(*) !<The number of non-zeros in each row
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATESEQAIJ",ERR,ERROR,*999)

    CALL MatCreateSeqAIJ(COMMUNICATOR,M,N,NUMBER_NZ_PERROW,NUMBER_NZ_EACHROW,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateSeqAIJ",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATESEQAIJ")
    RETURN
999 CALL ERRORS("PETSC_MATCREATESEQAIJ",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATESEQAIJ")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATESEQAIJ
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateSeqDense routine.
  SUBROUTINE PETSC_MATCREATESEQDENSE(COMMUNICATOR,M,N,MATRIX_DATA,A,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI Communicator
    INTEGER(INTG), INTENT(IN) :: M !<The number of rows
    INTEGER(INTG), INTENT(IN) :: N !<The number of columns
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*) !<Optional, the allocated matrix data
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATCREATESEQDENSE",ERR,ERROR,*999)

    CALL MatCreateSeqDense(COMMUNICATOR,M,N,MATRIX_DATA,A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatCreateSeqDense",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATCREATESEQDENSE")
    RETURN
999 CALL ERRORS("PETSC_MATCREATESEQDENSE",ERR,ERROR)
    CALL EXITS("PETSC_MATCREATESEQDENSE")
    RETURN 1
  END SUBROUTINE PETSC_MATCREATESEQDENSE

  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatSetType routine.
  SUBROUTINE PETSC_MATSETTYPE(A,MATRIXTYPE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix
    MatType, INTENT(IN) :: MATRIXTYPE !<The matrix type 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETTYPE",ERR,ERROR,*999)

    CALL MatSetType(A%MAT,MATRIXTYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETTYPE")
    RETURN
999 CALL ERRORS("PETSC_MATSETTYPE",ERR,ERROR)
    CALL EXITS("PETSC_MATSETTYPE")
    RETURN 1
  END SUBROUTINE PETSC_MATSETTYPE

  !
  !================================================================================================================================
  !    

  !>Buffer routine to the PETSc MatDestroy routine.
  SUBROUTINE PETSC_MATDESTROY(A,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATDESTROY",ERR,ERROR,*999)

    CALL MatDestroy(A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatDestroy",ERR,ERROR,*999)
    ENDIF
    A%MAT=PETSC_NULL_INTEGER
     
    CALL EXITS("PETSC_MATDESTROY")
    RETURN
999 CALL ERRORS("PETSC_MATDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_MATDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_MATDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringCreate routine.
  SUBROUTINE PETSC_MATFDCOLORINGCREATE(A,ISCOLORING,FDCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The PETSc matrix to create the FD coloring for
    TYPE(PETSC_ISCOLORING_TYPE), INTENT(IN) :: ISCOLORING !<The index set coloring to create the finite difference coloring for
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(OUT) :: FDCOLORING !<On exit, the matrix finite difference coloring
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGCREATE",ERR,ERROR,*999)

    CALL MatFDColoringCreate(A%MAT,ISCOLORING%ISCOLORING,FDCOLORING%MATFDCOLORING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatFDColoringCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATFDCOLORINGCREATE")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGCREATE",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGCREATE")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringDestroy routine.
  SUBROUTINE PETSC_MATFDCOLORINGDESTROY(MATFDCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The matrix finite difference coloring to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGDESTROY",ERR,ERROR,*999)

    CALL MatFDColoringDestroy(MATFDCOLORING%MATFDCOLORING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatFDColoringDestroy",ERR,ERROR,*999)
    ENDIF
    MATFDCOLORING%MATFDCOLORING=PETSC_NULL_INTEGER
     
    CALL EXITS("PETSC_MATFDCOLORINGDESTROY")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGDESTROY
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc MatFDColoring structure and destroy the MatFDColoring
  SUBROUTINE PETSC_MATFDCOLORINGFINALISE(MATFDCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The MatFDColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGFINALISE",ERR,ERROR,*999)

    IF(MATFDCOLORING%MATFDCOLORING/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_MATFDCOLORINGDESTROY(MATFDCOLORING,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATFDCOLORINGFINALISE")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc MatFDColoring structure
  SUBROUTINE PETSC_MATFDCOLORINGINITIALISE(MATFDCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The MatFDColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGINITIALISE",ERR,ERROR,*999)

    MATFDCOLORING%MATFDCOLORING=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_MATFDCOLORINGINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetFromOptions routine.
  SUBROUTINE PETSC_MATFDCOLORINGSETFROMOPTIONS(MATFDCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The matrix finite difference coloring to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGSETFROMOPTIONS",ERR,ERROR,*999)

    CALL MatFDColoringSetFromOptions(MATFDCOLORING%MATFDCOLORING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatFDColoringSetFromOptions",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATFDCOLORINGSETFROMOPTIONS")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGSETFROMOPTIONS",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGSETFROMOPTIONS")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetParameters routine.
  SUBROUTINE PETSC_MATFDCOLORINGSETPARAMETERS(MATFDCOLORING,RERROR,UMIN,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The matrix finite difference coloring to set
    REAL(DP) :: RERROR !<The square root of the relative error
    REAL(DP) :: UMIN !<MatFDColoring umin option
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGSETPARAMETERS",ERR,ERROR,*999)

    CALL MatFDColoringSetParameters(MATFDCOLORING%MATFDCOLORING,RERROR,UMIN,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatFDColoringSetParameters",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("PETSC_MATFDCOLORINGSETPARAMETERS")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGSETPARAMETERS",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGSETPARAMETERS")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGSETPARAMETERS

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetFunction routine.
  SUBROUTINE PETSC_MATFDCOLORINGSETFUNCTION(MATFDCOLORING,FFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The matrix finite difference coloring to set
    EXTERNAL FFUNCTION !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGSETFUNCTION",ERR,ERROR,*999)

!!Temporary until we move over to PETSc 3.
#if ( PETSC_VERSION_MAJOR == 3 )
    CALL MatFDColoringSetFunction(MATFDCOLORING%MATFDCOLORING,FFUNCTION,CTX,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatFDColoringSetFunction",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATFDCOLORINGSETFUNCTION")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGSETFUNCTION",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGSETFUNCTION")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGSETFUNCTION
        
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetFunctionSNES routine.
  SUBROUTINE PETSC_MATFDCOLORINGSETFUNCTIONSNES(MATFDCOLORING,FFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MATFDCOLORING_TYPE), INTENT(INOUT) :: MATFDCOLORING !<The matrix finite difference coloring to set
    EXTERNAL FFUNCTION !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFDCOLORINGSETFUNCTIONSNES",ERR,ERROR,*999)

#if ( PETSC_VERSION_MAJOR == 2 )
    CALL MatFDColoringSetFunctionSNES(MATFDCOLORING%MATFDCOLORING,FFUNCTION,CTX,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatFDColoringSetFunctionSNES",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATFDCOLORINGSETFUNCTIONSNES")
    RETURN
999 CALL ERRORS("PETSC_MATFDCOLORINGSETFUNCTIONSNES",ERR,ERROR)
    CALL EXITS("PETSC_MATFDCOLORINGSETFUNCTIONSNES")
    RETURN 1
  END SUBROUTINE PETSC_MATFDCOLORINGSETFUNCTIONSNES
        
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetArray routine.
  SUBROUTINE PETSC_MATGETARRAY(A,ARRAY,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT), TARGET :: A !<The matrix to get the array for
    REAL(DP), POINTER :: ARRAY(:) !<On exit, a pointer to the matrix array
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATGETARRAY",ERR,ERROR,*999)

    IF(ASSOCIATED(ARRAY)) THEN
      CALL FLAG_ERROR("Array is already associated",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("PETSC_MATGETARRAY not implemented",ERR,ERROR,*999)
      !CALL MatGetArray(A%MAT,A%MAT_DATA,A%MAT_OFFSET,ERR)
      !IF(ERR/=0) THEN
      !  IF(PETSC_HANDLE_ERROR) THEN
      !    CHKERRQ(ERR)
      !  ENDIF
      !  CALL FLAG_ERROR("PETSc error in MatGetArray",ERR,ERROR,*999)
      !ENDIF
    ENDIF
    
    CALL EXITS("PETSC_MATGETARRAY")
    RETURN
999 CALL ERRORS("PETSC_MATGETARRAY",ERR,ERROR)
    CALL EXITS("PETSC_MATGETARRAY")
    RETURN 1
  END SUBROUTINE PETSC_MATGETARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetArrayF90 routine.
  SUBROUTINE PETSC_MATGETARRAYF90(A,ARRAY,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT), TARGET :: A !<The matrix to get the array for
    REAL(DP), POINTER :: ARRAY(:,:) !<On exit, a pointer to the matrix array
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATGETARRAYF90",ERR,ERROR,*999)

    IF(ASSOCIATED(ARRAY)) THEN
      CALL FLAG_ERROR("Array is already associated",ERR,ERROR,*999)
    ELSE
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 3 )
      CALL MatSeqAIJGetArrayF90(A%MAT,ARRAY,ERR)
#else
      CALL MatGetArrayF90(A%MAT,ARRAY,ERR)
#endif
      IF(ERR/=0) THEN
        IF(PETSC_HANDLE_ERROR) THEN
          CHKERRQ(ERR)
        ENDIF
        CALL FLAG_ERROR("PETSc error in MatGetArrayF90",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    CALL EXITS("PETSC_MATGETARRAYF90")
    RETURN
999 CALL ERRORS("PETSC_MATGETARRAYF90",ERR,ERROR)
    CALL EXITS("PETSC_MATGETARRAYF90")
    RETURN 1
  END SUBROUTINE PETSC_MATGETARRAYF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetColoring routine.
  SUBROUTINE PETSC_MATGETCOLORING(A,COLORING_TYPE,ISCOLORING,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to get the ownership range of
    MatColoringType, INTENT(IN) :: COLORING_TYPE !<The coloring type
    TYPE(PETSC_ISCOLORING_TYPE), INTENT(OUT) :: ISCOLORING !<On exit, the index set coloring
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    MatColoring :: COLORING
#endif

    CALL ENTERS("PETSC_MATGETCOLORING",ERR,ERROR,*999)
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    CALL MatColoringCreate(A%MAT,COLORING,ERR);
    CALL MatColoringSetType(COLORING,COLORING_TYPE,ERR);
    CALL MatColoringSetFromOptions(COLORING,ERR);
    CALL MatColoringApply(COLORING,ISCOLORING%ISCOLORING,ERR);
    CALL MatColoringDestroy(COLORING,ERR);
#else
    CALL MatGetColoring(A%MAT,COLORING_TYPE,ISCOLORING%ISCOLORING,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatGetColoring",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATGETCOLORING")
    RETURN
999 CALL ERRORS("PETSC_MATGETCOLORING",ERR,ERROR)
    CALL EXITS("PETSC_MATGETCOLORING")
    RETURN 1
  END SUBROUTINE PETSC_MATGETCOLORING

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetInfo routine.
  SUBROUTINE PETSC_MATGETINFO(A,FLAG,INFO,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to get the information of
    MatInfoType, INTENT(IN) :: FLAG !<The matrix information collective to return
    MatInfo, INTENT(OUT) :: INFO(MAT_INFO_SIZE) !<On return, the matrix information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATGETINFO",ERR,ERROR,*999)

    CALL MatGetInfo(A%MAT,FLAG,INFO,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatGetInfo",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATGETINFO")
    RETURN
999 CALL ERRORS("PETSC_MATGETINFO",ERR,ERROR)
    CALL EXITS("PETSC_MATGETINFO")
    RETURN 1
  END SUBROUTINE PETSC_MATGETINFO
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetOwnershipRange routine.
  SUBROUTINE PETSC_MATGETOWNERSHIPRANGE(A,FIRST_ROW,LAST_ROW,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to get the ownership range of
    INTEGER(INTG), INTENT(OUT) :: FIRST_ROW !<On exit, the first row for the matrix
    INTEGER(INTG), INTENT(OUT) :: LAST_ROW !<On exit, the last row for the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATGETOWNERSHIPRANGE",ERR,ERROR,*999)

    CALL MatGetOwnershipRange(A%MAT,FIRST_ROW,LAST_ROW,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatGetOwnershipRange",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATGETOWNERSHIPRANGE")
    RETURN
999 CALL ERRORS("PETSC_MATGETOWNERSHIPRANGE",ERR,ERROR)
    CALL EXITS("PETSC_MATGETOWNERSHIPRANGE")
    RETURN 1
  END SUBROUTINE PETSC_MATGETOWNERSHIPRANGE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetRow routine.
  SUBROUTINE PETSC_MATGETROW(A,ROW_NUMBER,NUMBER_OF_COLUMNS,COLUMNS,VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to get the array for
    INTEGER(INTG), INTENT(IN) :: ROW_NUMBER !<The row number to get the row values for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COLUMNS !<On return, the number of nonzero columns in the row
    INTEGER(INTG), INTENT(OUT) :: COLUMNS(:) !<On return, the column numbers for the nonzero columns in the row
    REAL(DP), INTENT(OUT) :: VALUES(:) !<On exit, the nonzero values in the row.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATGETROW",ERR,ERROR,*999)

    CALL MatGetRow(A%MAT,ROW_NUMBER,NUMBER_OF_COLUMNS,COLUMNS,VALUES,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatGetRow",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATGETROW")
    RETURN
999 CALL ERRORS("PETSC_MATGETROW",ERR,ERROR)
    CALL EXITS("PETSC_MATGETROW")
    RETURN 1
  END SUBROUTINE PETSC_MATGETROW
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetValues routine.
  SUBROUTINE PETSC_MATGETVALUES(A,M,M_INDICES,N,N_INDICES,VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to get the values of
    INTEGER(INTG), INTENT(IN) :: M !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: N !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*) !<The column indices
    REAL(DP), INTENT(OUT) :: VALUES(*) !<The values to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATGETVALUES",ERR,ERROR,*999)

    CALL MatGetValues(A%MAT,M,M_INDICES,N,N_INDICES,VALUES,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatGetValues",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATGETVALUES")
    RETURN
999 CALL ERRORS("PETSC_MATGETVALUES",ERR,ERROR)
    CALL EXITS("PETSC_MATGETVALUES")
    RETURN 1
  END SUBROUTINE PETSC_MATGETVALUES
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc Mat structure and destroy the KSP
  SUBROUTINE PETSC_MATFINALISE(MAT_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: MAT_ !<The MAT to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATFINALISE",ERR,ERROR,*999)

    IF(MAT_%MAT/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_MATDESTROY(MAT_,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATFINALISE")
    RETURN
999 CALL ERRORS("PETSC_MATFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_MATFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_MATFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Mat structure
  SUBROUTINE PETSC_MATINITIALISE(MAT_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: MAT_ !<The MAT to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATINITIALISE",ERR,ERROR,*999)

    MAT_%MAT=PETSC_NULL_INTEGER
    !MAT_%MAT_DATA(1)=0
    !MAT_%MAT_OFFSET=0
    
    CALL EXITS("PETSC_MATINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_MATINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_MATINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_MATINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatRestoreArray routine.
  SUBROUTINE PETSC_MATRESTOREARRAY(A,ARRAY,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to restore the array for
    REAL(DP), POINTER :: ARRAY(:) !<A pointer to the matrix array to restore
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATRESTOREARRAY",ERR,ERROR,*999)

    CALL FLAG_ERROR("PETSC_MATRESTOREARRAY not implemented",ERR,ERROR,*999)
    !CALL MatRestoreArray(A%MAT,A%MAT_DATA,A%MAT_OFFSET,ERR)
    !IF(ERR/=0) THEN
    !  IF(PETSC_HANDLE_ERROR) THEN
    !    CHKERRQ(ERR)
    !  ENDIF
    !  CALL FLAG_ERROR("PETSc error in MatRestoreArray",ERR,ERROR,*999)
    !ENDIF
    
    CALL EXITS("PETSC_MATRESTOREARRAY")
    RETURN
999 CALL ERRORS("PETSC_MATRESTOREARRAY",ERR,ERROR)
    CALL EXITS("PETSC_MATRESTOREARRAY")
    RETURN 1
  END SUBROUTINE PETSC_MATRESTOREARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatRestoreArrayF90 routine.
  SUBROUTINE PETSC_MATRESTOREARRAYF90(A,ARRAY,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to restore the array for
    REAL(DP), POINTER :: ARRAY(:,:) !<A pointer to the matrix array to restore
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATRESTOREARRAYF90",ERR,ERROR,*999)

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 3 )
    CALL MatSeqAIJRestoreArrayF90(A%MAT,ARRAY,ERR)
#else
    CALL MatRestoreArrayF90(A%MAT,ARRAY,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatRestoreArrayF90",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATRESTOREARRAYF90")
    RETURN
999 CALL ERRORS("PETSC_MATRESTOREARRAYF90",ERR,ERROR)
    CALL EXITS("PETSC_MATRESTOREARRAYF90")
    RETURN 1
  END SUBROUTINE PETSC_MATRESTOREARRAYF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatRestoreRow routine.
  SUBROUTINE PETSC_MATRESTOREROW(A,ROW_NUMBER,NUMBER_OF_COLUMNS,COLUMNS,VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to restore the row for
    INTEGER(INTG), INTENT(IN) :: ROW_NUMBER !<The row number to restore the row values for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COLUMNS !<The number of nonzero columns in the row
    INTEGER(INTG), INTENT(OUT) :: COLUMNS(:) !<The column numbers for the nonzero columns in the row
    REAL(DP), INTENT(OUT) :: VALUES(:) !<The nonzero values in the row to restore.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATRESTOREROW",ERR,ERROR,*999)

    CALL MatRestoreRow(A%MAT,ROW_NUMBER,NUMBER_OF_COLUMNS,COLUMNS,VALUES,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatRestoreRow.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATRESTOREROW")
    RETURN
999 CALL ERRORS("PETSC_MATRESTOREROW",ERR,ERROR)
    CALL EXITS("PETSC_MATRESTOREROW")
    RETURN 1
  END SUBROUTINE PETSC_MATRESTOREROW
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetLocalToGlobalMapping routine.
  SUBROUTINE PETSC_MATSETLOCALTOGLOBALMAPPING(A,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the local to global mapping for
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(IN) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETLOCALTOGLOBALMAPPING",ERR,ERROR,*999)

    CALL MatSetLocalToGlobalMapping(A%MAT,CTX%ISLOCALTOGLOBALMAPPING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetLocalToGlobalMapping",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETLOCALTOGLOBALMAPPING")
    RETURN
999 CALL ERRORS("PETSC_MATSETLOCALTOGLOBALMAPPING",ERR,ERROR)
    CALL EXITS("PETSC_MATSETLOCALTOGLOBALMAPPING")
    RETURN 1
  END SUBROUTINE PETSC_MATSETLOCALTOGLOBALMAPPING
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetOption routine.
  SUBROUTINE PETSC_MATSETOPTION(A,OPTION,FLAG,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the option for
    MatOption, INTENT(IN) :: OPTION !<The option to set
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
    PetscBool, INTENT(IN) :: FLAG !<The option flag to set
#else
    PetscTruth, INTENT(IN) :: FLAG !<The option flag to set
#endif    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETOPTION",ERR,ERROR,*999)

    CALL MatSetOption(A%MAT,OPTION,FLAG,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetOption",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETOPTION")
    RETURN
999 CALL ERRORS("PETSC_MATSETOPTION",ERR,ERROR)
    CALL EXITS("PETSC_MATSETOPTION")
    RETURN 1
  END SUBROUTINE PETSC_MATSETOPTION
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetSizes routine.
  SUBROUTINE PETSC_MATSETSIZES(A,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the size of
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<Number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<Number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<Number of global rows
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<Number of global columns
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETSIZES",ERR,ERROR,*999)

    CALL MatSetSizes(A%MAT,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetSizes",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETSIZES")
    RETURN
999 CALL ERRORS("PETSC_MATSETSIZES",ERR,ERROR)
    CALL EXITS("PETSC_MATSETSIZES")
    RETURN 1
  END SUBROUTINE PETSC_MATSETSIZES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValue routine.
  SUBROUTINE PETSC_MATSETVALUE(A,ROW,COL,VALUE,INSERT_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: ROW !<The row index
    INTEGER(INTG), INTENT(IN) :: COL !<The column index
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETVALUE",ERR,ERROR,*999)

    CALL MatSetValue(A%MAT,ROW,COL,VALUE,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetValue",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETVALUE")
    RETURN
999 CALL ERRORS("PETSC_MATSETVALUE",ERR,ERROR)
    CALL EXITS("PETSC_MATSETVALUE")
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValues routine.
  SUBROUTINE PETSC_MATSETVALUES(A,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: M !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: N !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*) !<The column indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETVALUES",ERR,ERROR,*999)

    CALL MatSetValues(A%MAT,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetValues",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETVALUES")
    RETURN
999 CALL ERRORS("PETSC_MATSETVALUES",ERR,ERROR)
    CALL EXITS("PETSC_MATSETVALUES")
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValueLocal routine.
  SUBROUTINE PETSC_MATSETVALUELOCAL(A,ROW,COL,VALUE,INSERT_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: ROW !<The row index
    INTEGER(INTG), INTENT(IN) :: COL !<The column index
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETVALUELOCAL",ERR,ERROR,*999)

    CALL MatSetValueLocal(A%MAT,ROW,COL,VALUE,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetValueLocal",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETVALUELOCAL")
    RETURN
999 CALL ERRORS("PETSC_MATSETVALUELOCAL",ERR,ERROR)
    CALL EXITS("PETSC_MATSETVALUELOCAL")
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUELOCAL
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValuesLocal routine.
  SUBROUTINE PETSC_MATSETVALUESLOCAL(A,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: M !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: N !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*) !<The column indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATSETVALUESLOCAL",ERR,ERROR,*999)

    CALL MatSetValuesLocal(A%MAT,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatSetValuesLocal",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATSETVALUESLOCAL")
    RETURN
999 CALL ERRORS("PETSC_MATSETVALUESLOCAL",ERR,ERROR)
    CALL EXITS("PETSC_MATSETVALUESLOCAL")
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUESLOCAL
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatView routine.
  SUBROUTINE PETSC_MATVIEW(A,V,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to view
    PetscViewer, INTENT(IN) :: V !<The viewer
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATVIEW",ERR,ERROR,*999)

    CALL MatView(A%MAT,V,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatView",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATVIEW")
    RETURN
999 CALL ERRORS("PETSC_MATVIEW",ERR,ERROR)
    CALL EXITS("PETSC_MATVIEW")
    RETURN 1
  END SUBROUTINE PETSC_MATVIEW
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatZeroEntries routine.
  SUBROUTINE PETSC_MATZEROENTRIES(A,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The matrix to zero the entries of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_MATZEROENTRIES",ERR,ERROR,*999)

    CALL MatZeroEntries(A%MAT,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatZeroEntries",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_MATZEROENTRIES")
    RETURN
999 CALL ERRORS("PETSC_MATZEROENTRIES",ERR,ERROR)
    CALL EXITS("PETSC_MATZEROENTRIES")
    RETURN 1
  END SUBROUTINE PETSC_MATZEROENTRIES
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc PC structure and destroy the KSP
  SUBROUTINE PETSC_PCFINALISE(PC_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<The PC to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_PCFINALISE",ERR,ERROR,*999)

    IF(PC_%PC_/=PETSC_NULL_INTEGER) THEN
      !Do nothing - should be destroyed when the KSP is destroyed.
    ENDIF
    
    CALL EXITS("PETSC_PCFINALISE")
    RETURN
999 CALL ERRORS("PETSC_PCFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_PCFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_PCFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc PC structure
  SUBROUTINE PETSC_PCINITIALISE(PC_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<The PC to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_PCINITIALISE",ERR,ERROR,*999)

    PC_%PC_=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_PCINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_PCINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_PCINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_PCINITIALISE
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR == 3 )
  !>Buffer routine to the PETSc PCFactoSetMatSolverPackage routine.
  SUBROUTINE PETSC_PCFACTORSETMATSOLVERPACKAGE(PC_,SOLVER_PACKAGE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<The preconditioner to set the solver package for
    MatSolverPackage, INTENT(IN) :: SOLVER_PACKAGE !<The solver package to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_PCFACTORSETMATSOLVERPACKAGE",ERR,ERROR,*999)

    CALL PCFactorSetMatSolverPackage(PC_%PC_,SOLVER_PACKAGE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PCFactorSetMatSolverPackage",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_PCFACTORSETMATSOLVERPACKAGE")
    RETURN
999 CALL ERRORS("PETSC_PCFACTORSETMATSOLVERPACKAGE",ERR,ERROR)
    CALL EXITS("PETSC_PCFACTORSETMATSOLVERPACKAGE")
    RETURN 1
  END SUBROUTINE PETSC_PCFACTORSETMATSOLVERPACKAGE
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR == 3 )
  !>Buffer routine to the PETSc PCFactorSetUpMatSolverPackage routine.
  SUBROUTINE Petsc_PCFactorSetUpMatSolverPackage(PC_,err,error,*)

    !Argument Variables
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<The preconditioner to set the solver package for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("Petsc_PCFactorSetUpMatSolverPackage",err,error,*999)

    CALL PCFactorSetUpMatSolverPackage(PC_%PC_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PCFactorSetUpMatSolverPackage",err,error,*999)
    ENDIF
    
    CALL EXITS("Petsc_PCFactorSetUpMatSolverPackage")
    RETURN
999 CALL ERRORS("Petsc_PCFactorSetUpMatSolverPackage",err,error)
    CALL EXITS("Petsc_PCFactorSetUpMatSolverPackage")
    RETURN 1
  END SUBROUTINE Petsc_PCFactorSetUpMatSolverPackage
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR == 3 )
  !>Buffer routine to the PETSc PCFactorGetMatrix routine.
  SUBROUTINE Petsc_PCFactorGetMatrix(PC_,factoredMatrix,err,error,*)

    !Argument Variables
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<The preconditioner to set the solver package for
    TYPE(PETSC_MAT_TYPE), INTENT(OUT) :: factoredMatrix !<The factored matrix to get from preconditioner context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("Petsc_PCFactorGetMatrix",err,error,*999)

    CALL PCFactorGetMatrix(PC_%PC_,factoredMatrix,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PCFactorGetMatrix",err,error,*999)
    ENDIF
    
    CALL EXITS("Petsc_PCFactorGetMatrix")
    RETURN
999 CALL ERRORS("Petsc_PCFactorGetMatrix",err,error)
    CALL EXITS("Petsc_PCFactorGetMatrix")
    RETURN 1
  END SUBROUTINE Petsc_PCFactorGetMatrix
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 )
  !>Buffer routine to the PETSc MatMumpsSetIcntl routine.
  SUBROUTINE Petsc_MatMumpsSetIcntl(factoredMatrix,icntl,ival,err,error,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: factoredMatrix !<The factored matrix from PETSc-MUMPS interface
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS ICNTL integer control parameter
    INTEGER(INTG), INTENT(IN) :: ival !<The MUMPS ICNTL integer value to set: ICNTL(icntl)=ival
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("Petsc_MatMumpsSetIcntl",err,error,*999)

    CALL MatMumpsSetIcntl(factoredMatrix,icntl,ival,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatMumpsSetIcntl",err,error,*999)
    ENDIF
    
    CALL EXITS("Petsc_MatMumpsSetIcntl")
    RETURN
999 CALL ERRORS("Petsc_MatMumpsSetIcntl",err,error)
    CALL EXITS("Petsc_MatMumpsSetIcntl")
    RETURN 1
  END SUBROUTINE Petsc_MatMumpsSetIcntl
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 4 )
  !>Buffer routine to the PETSc MatMumpsSetCntl routine.
  SUBROUTINE Petsc_MatMumpsSetCntl(factoredMatrix,icntl,val,err,error,*)

    !Argument Variables
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: factoredMatrix !<The factored matrix from PETSc-MUMPS interface
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS CNTL integer control parameter
    REAL(DP), INTENT(IN) :: val !<The MUMPS CNTL real value to set: CNTL(icntl)=val
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("Petsc_MatMumpsSetCntl",err,error,*999)

    CALL MatMumpsSetCntl(factoredMatrix,icntl,val,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in MatMumpsSetCntl",err,error,*999)
    ENDIF
    
    CALL EXITS("Petsc_MatMumpsSetCntl")
    RETURN
999 CALL ERRORS("Petsc_MatMumpsSetCntl",err,error)
    CALL EXITS("Petsc_MatMumpsSetCntl")
    RETURN 1
  END SUBROUTINE Petsc_MatMumpsSetCntl
#endif
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCSetType routine.
  SUBROUTINE PETSC_PCSETTYPE(PC_,METHOD,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_PC_TYPE), INTENT(INOUT) :: PC_ !<The preconditioner to set the type of
    PCType, INTENT(IN) :: METHOD !<The preconditioning method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_PCSETTYPE",ERR,ERROR,*999)

    CALL PCSetType(PC_%PC_,METHOD,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in PCSetType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_PCSETTYPE")
    RETURN
999 CALL ERRORS("PETSC_PCSETTYPE",ERR,ERROR)
    CALL EXITS("PETSC_PCSETTYPE")
    RETURN 1
  END SUBROUTINE PETSC_PCSETTYPE
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc SNES structure and destroy the SNES
  SUBROUTINE PETSC_SNESFINALISE(SNES_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESFINALISE",ERR,ERROR,*999)

    IF(SNES_%SNES_/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_SNESDESTROY(SNES_,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESFINALISE")
    RETURN
999 CALL ERRORS("PETSC_SNESFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_SNESFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_SNESFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc SNES structure
  SUBROUTINE PETSC_SNESINITIALISE(SNES_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The snes to 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESINITIALISE",ERR,ERROR,*999)

    SNES_%SNES_=PETSC_NULL_INTEGER
     
    CALL EXITS("PETSC_SNESINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_SNESINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_SNESINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_SNESINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESCreate routine.
  SUBROUTINE PETSC_SNESCREATE(COMMUNICATOR,SNES_,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator for the SNES creation
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<On exit, the SNES information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESCREATE",ERR,ERROR,*999)

    CALL SNESCreate(COMMUNICATOR,SNES_%SNES_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESCREATE")
    RETURN
999 CALL ERRORS("PETSC_SNESCREATE",ERR,ERROR)
    CALL EXITS("PETSC_SNESCREATE")
    RETURN 1
  END SUBROUTINE PETSC_SNESCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESDestroy routine.
  SUBROUTINE PETSC_SNESDESTROY(SNES_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESDESTROY",ERR,ERROR,*999)

    CALL SNESDestroy(SNES_%SNES_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESDestroy",ERR,ERROR,*999)
    ENDIF
    SNES_%SNES_=PETSC_NULL_INTEGER
     
    CALL EXITS("PETSC_SNESDESTROY")
    RETURN
999 CALL ERRORS("PETSC_SNESDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_SNESDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_SNESDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetConvergedReason routine.
  SUBROUTINE PETSC_SNESGETCONVERGEDREASON(SNES_,REASON,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to get the converged reason for
    INTEGER(INTG), INTENT(OUT) :: REASON !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESGETCONVERGEDREASON",ERR,ERROR,*999)

    CALL SNESGetConvergedReason(SNES_%SNES_,REASON,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetConvergedReason",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESGETCONVERGEDREASON")
    RETURN
999 CALL ERRORS("PETSC_SNESGETCONVERGEDREASON",ERR,ERROR)
    CALL EXITS("PETSC_SNESGETCONVERGEDREASON")
    RETURN 1
  END SUBROUTINE PETSC_SNESGETCONVERGEDREASON
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetFunctionNorm routine.
  SUBROUTINE PETSC_SNESGETFUNCTIONNORM(SNES_,FUNCTION_NORM,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to get the function norm for
    REAL(DP), INTENT(OUT) :: FUNCTION_NORM !<On exit, the function norm
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESGETFUNCTIONNORM",ERR,ERROR,*999)
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    CALL SNESGetFunction(SNES_%SNES_,FUNCTION_NORM,ERR)
    CALL VecNorm(SNES_%SNES_,FUNCTION_NORM,ERR)
#else
    CALL SNESGetFunctionNorm(SNES_%SNES_,FUNCTION_NORM,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetFunctionNorm",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESGETFUNCTIONNORM")
    RETURN
999 CALL ERRORS("PETSC_SNESGETFUNCTIONNORM",ERR,ERROR)
    CALL EXITS("PETSC_SNESGETFUNCTIONNORM")
    RETURN 1
  END SUBROUTINE PETSC_SNESGETFUNCTIONNORM

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSC SNESGetSolutionUpdate routine.
  SUBROUTINE Petsc_SnesGetSolutionUpdate(snes_,solutionUpdate,err,error,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: snes_ !<The SNES to get the solution update for
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: solutionUpdate !<On exit, the solution update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("Petsc_SnesGetSolutionUpdate",err,error,*999)

    CALL SNESGetSolutionUpdate(snes_%SNES_,solutionUpdate%VEC,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetSolutionUpdate",err,error,*999)
    ENDIF

    CALL EXITS("Petsc_SnesGetSolutionUpdate")
    RETURN
999 CALL ERRORS("Petsc_SnesGetSolutionUpdate",err,error)
    CALL EXITS("Petsc_SnesGetSolutionUpdate")
    RETURN 1
  END SUBROUTINE Petsc_SnesGetSolutionUpdate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunctionNorm routine.
  SUBROUTINE PETSC_SNESSETFUNCTIONNORM(SNES_,FUNCTION_NORM,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to get the function norm for
    REAL(DP), INTENT(OUT) :: FUNCTION_NORM !<On exit, the function norm
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETFUNCTIONNORM",ERR,ERROR,*999)
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    CALL SNESGetFunction(SNES_%SNES_,FUNCTION_NORM,ERR)
    CALL VecNorm(SNES_%SNES_,FUNCTION_NORM,ERR)
#else
    CALL SNESSetFunctionNorm(SNES_%SNES_,FUNCTION_NORM,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetFunctionNorm",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETFUNCTIONNORM")
    RETURN
999 CALL ERRORS("PETSC_SNESSETFUNCTIONNORM",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETFUNCTIONNORM")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETFUNCTIONNORM
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the petsc SnesLineSearchSetNorms routine.
  SUBROUTINE PETSC_SnesLineSearchSetNorms(SNES_,XNORM,FNORM,YNORM,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to get the computed norms for X, Y, and F
    REAL(DP), INTENT(INOUT) :: XNORM !<On exit, the norm of the current solution
    REAL(DP), INTENT(INOUT) :: FNORM !<On exit, the norm of the current function
    REAL(DP), INTENT(INOUT) :: YNORM !<On exit, the norm of the current update
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("petsc_SnesLineSearchSetNorms",ERR,ERROR,*999)

    CALL SnesLineSearchSetNorms(SNES_%SNES_,XNORM,FNORM,YNORM,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("petsc error in SnesLineSearchSetNorms",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("petsc_SnesLineSearchSetNorms")
    RETURN
999 CALL ERRORS("petsc_SnesLineSearchSetNorms",ERR,ERROR)
    CALL EXITS("petsc_SnesLineSearchSetNorms")
    RETURN 1
  END SUBROUTINE petsc_SnesLineSearchSetNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the petsc SnesLineSearchGetNorms routine.
  SUBROUTINE petsc_SnesLineSearchGetNorms(lineSearch,XNORM,FNORM,YNORM,ERR,ERROR,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to get the norms for X, Y, and F from.
    REAL(DP), INTENT(INOUT) :: XNORM !<On exit, the norm of the current solution
    REAL(DP), INTENT(INOUT) :: FNORM !<On exit, the norm of the current function
    REAL(DP), INTENT(INOUT) :: YNORM !<On exit, the norm of the current update
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("petsc_SnesLineSearchGetNorms",ERR,ERROR,*999)

    CALL SnesLineSearchGetNorms(lineSearch%snesLineSearch,XNORM,FNORM,YNORM,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("petsc error in SnesLineSearchGetNorms",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("petsc_SnesLineSearchGetNorms")
    RETURN
999 CALL ERRORS("petsc_SnesLineSearchGetNorms",ERR,ERROR)
    CALL EXITS("petsc_SnesLineSearchGetNorms")
    RETURN 1
  END SUBROUTINE petsc_SnesLineSearchGetNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetIterationNumber routine.
  SUBROUTINE PETSC_SNESGETITERATIONNUMBER(SNES_,ITERATION_NUMBER,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to get the iteration number for
    INTEGER(INTG), INTENT(OUT) :: ITERATION_NUMBER !<On exit, the number of iterations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESGETITERATIONNUMBER",ERR,ERROR,*999)

    CALL SNESGetIterationNumber(SNES_%SNES_,ITERATION_NUMBER,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetIterationNumber",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESGETITERATIONNUMBER")
    RETURN
999 CALL ERRORS("PETSC_SNESGETITERATIONNUMBER",ERR,ERROR)
    CALL EXITS("PETSC_SNESGETITERATIONNUMBER")
    RETURN 1
  END SUBROUTINE PETSC_SNESGETITERATIONNUMBER
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetKSP routine.
  SUBROUTINE PETSC_SNESGETKSP(SNES_,KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to get the iteration number for
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<On exit, the KSP to associated with the SNES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESGETKSP",ERR,ERROR,*999)

    CALL SNESGetKSP(SNES_%SNES_,KSP_%KSP_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetKSP",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESGETKSP")
    RETURN
999 CALL ERRORS("PETSC_SNESGETKSP",ERR,ERROR)
    CALL EXITS("PETSC_SNESGETKSP")
    RETURN 1
  END SUBROUTINE PETSC_SNESGETKSP

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !>Buffer routine to the PETSc SNESGetSNESLineSearch routine.
  SUBROUTINE Petsc_SnesGetSnesLineSearch(snes_,lineSearch,err,error,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: snes_ !<The SNES to get the SNES line search for
    TYPE(PetscSnesLinesearchType), INTENT(OUT) :: lineSearch !<On return, the SNES line search object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("Petsc_SnesGetSnesLineSearch",err,error,*999)
    
#if ( PETSC_VERSION_MINOR >= 4 )
    CALL SNESGetLineSearch(snes_%snes_,lineSearch%snesLineSearch,err)
#else
    CALL SNESGetSNESLineSearch(snes_%snes_,lineSearch%snesLineSearch,err)
#endif
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetLineSearch",err,error,*999)
    ENDIF

    CALL Exits("Petsc_SnesGetSnesLineSearch")
    RETURN
999 CALL Errors("Petsc_SnesGetSnesLineSearch",err,error)
    CALL Exits("Petsc_SnesGetSnesLineSearch")
    RETURN 1
  END SUBROUTINE Petsc_SnesGetSnesLineSearch
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !Finalise the PETSc SNES line search structure
  SUBROUTINE Petsc_SnesLineSearchFinalise(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("Petsc_SnesLineSearchFinalise",err,error,*999)

    ! We don't actually call PETSc's SNESLineSearchDestroy as PETSc accesses
    ! the LineSearch when calling SNESDestroy and also destroys it when
    ! calling SNESDestroy, so we'll just let PETSc clean it up.
    lineSearch%snesLineSearch=PETSC_NULL_INTEGER

    CALL Exits("Petsc_SnesLineSearchFinalise")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchFinalise",err,error)
    CALL Exits("Petsc_SnesLineSearchFinalise")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchFinalise

  !
  !================================================================================================================================
  !

  !Initialise the PETSc SNES line search structure
  SUBROUTINE Petsc_SnesLineSearchInitialise(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("Petsc_SnesLineSearchInitialise",err,error,*999)

    lineSearch%snesLineSearch=PETSC_NULL_INTEGER

    CALL Exits("Petsc_SnesLineSearchInitialise")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchInitialise",err,error)
    CALL Exits("Petsc_SnesLineSearchInitialise")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchInitialise
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
  !>Buffer routine to the PETSc SNESLineSearchSet routine.
  SUBROUTINE PETSC_SNESLINESEARCHSET(SNES_,LINESEARCH_TYPE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the line search for
    INTEGER(INTG), INTENT(IN) :: LINESEARCH_TYPE !<The line search type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESLINESEARCHSET",ERR,ERROR,*999)

    SELECT CASE(LINESEARCH_TYPE)
    CASE(PETSC_SNES_LINESEARCH_NONORMS)
      CALL SNESLineSearchSet(SNES_%SNES_,SNESLINESEARCHNONORMS,PETSC_NULL_OBJECT,ERR)
    CASE(PETSC_SNES_LINESEARCH_NO)
      CALL SNESLineSearchSet(SNES_%SNES_,SNESLINESEARCHNO,PETSC_NULL_OBJECT,ERR)
    CASE(PETSC_SNES_LINESEARCH_QUADRATIC)
      CALL SNESLineSearchSet(SNES_%SNES_,SNESLINESEARCHQUADRATIC,PETSC_NULL_OBJECT,ERR)
    CASE(PETSC_SNES_LINESEARCH_CUBIC)
      CALL SNESLineSearchSet(SNES_%SNES_,SNESLINESEARCHCUBIC,PETSC_NULL_OBJECT,ERR)      
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid line search type.",ERR,ERROR,*999)
    END SELECT
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESLineSearchSet.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESLINESEARCHSET")
    RETURN
999 CALL ERRORS("PETSC_SNESLINESEARCHSET",ERR,ERROR)
    CALL EXITS("PETSC_SNESLINESEARCHSET")
    RETURN 1
  END SUBROUTINE PETSC_SNESLINESEARCHSET
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2 )
  !>Buffer routine to the PETSc SNESLineSearchSetMonitor routine.
  SUBROUTINE Petsc_SnesLineSearchSetMonitor(lineSearch,monitorLinesearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set whether to output linesearch debug information
    PetscBool, INTENT(IN) :: monitorLinesearch !<Whether to output linesearch debug information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("Petsc_SnesLineSearchSetMonitor",err,error,*999)

    CALL SnesLineSearchSetMonitor(lineSearch%snesLineSearch,monitorLinesearch,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetMonitor",err,error,*999)
    ENDIF

    CALL Exits("Petsc_SnesLineSearchSetMonitor")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchSetMonitor",err,error)
    CALL Exits("Petsc_SnesLineSearchSetMonitor")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetMonitor
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !>Buffer routine to the PETSc SNESLineSearchSetComputeNorms routine.
  SUBROUTINE Petsc_SnesLineSearchSetComputeNorms(lineSearch,computeNorms,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set whether to compute norms
    PetscBool, INTENT(IN) :: computeNorms !<Whether to compute norms
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("Petsc_SnesLineSearchSetComputeNorms",err,error,*999)

    CALL SNESLineSearchSetComputeNorms(lineSearch%snesLineSearch,computeNorms,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetComputeNorms",err,error,*999)
    ENDIF

    CALL Exits("Petsc_SnesLineSearchSetComputeNorms")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchSetComputeNorms",err,error)
    CALL Exits("Petsc_SnesLineSearchSetComputeNorms")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetComputeNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchComputeNorms routine.
  SUBROUTINE Petsc_SnesLineSearchComputeNorms(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to compute norms for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("Petsc_SnesLineSearchComputeNorms",err,error,*999)

    CALL SnesLineSearchComputeNorms(lineSearch%snesLineSearch,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SnesLineSearchComputeNorms",err,error,*999)
    ENDIF

    CALL Exits("Petsc_SnesLineSearchComputeNorms")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchComputeNorms",err,error)
    CALL Exits("Petsc_SnesLineSearchComputeNorms")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchComputeNorms
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !>Buffer routine to the PETSc SNESLineSearchSetOrder routine.
  SUBROUTINE Petsc_SnesLineSearchSetOrder(lineSearch,lineSearchOrder,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the line search order for
    SNESLineSearchOrder, INTENT(IN) :: lineSearchOrder !<The line search order to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("Petsc_SnesLineSearchSetOrder",err,error,*999)

    SELECT CASE(LINESEARCHORDER)
    CASE(PETSC_SNES_LINESEARCH_LINEAR)
      CALL SNESLineSearchSetOrder(lineSearch%snesLineSearch,PETSC_SNES_LINESEARCH_LINEAR,err)
    CASE(PETSC_SNES_LINESEARCH_QUADRATIC)
      CALL SNESLineSearchSetOrder(lineSearch%snesLineSearch,PETSC_SNES_LINESEARCH_QUADRATIC,err)
    CASE(PETSC_SNES_LINESEARCH_CUBIC)
      CALL SNESLineSearchSetOrder(lineSearch%snesLineSearch,PETSC_SNES_LINESEARCH_CUBIC,err)
    CASE DEFAULT
      CALL FlagError("Invalid line search order.",err,error,*999)
    END SELECT
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetOrder",err,error,*999)
    ENDIF

    CALL Exits("Petsc_SnesLineSearchSetOrder")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchSetOrder",err,error)
    CALL Exits("Petsc_SnesLineSearchSetOrder")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetOrder
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 3 )
  !>Buffer routine to the PETSc SNESLineSearchSetParams routine.
  SUBROUTINE PETSC_SNESLINESEARCHSETPARAMS(SNES_,ALPHA,MAXSTEP,STEPTOL,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the line search parameters for
    REAL(DP), INTENT(IN) :: ALPHA !<The scalar such that 0.5f_{n+1} . f_{n+1} <= .5*f_n . f_n - alpha |f_n . J . f_n| 
    REAL(DP), INTENT(IN) :: MAXSTEP !<The maximum norm of the update vector
    REAL(DP), INTENT(IN) :: STEPTOL !<the minimum norm fraction of the the original step after scaling
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESLINESEARCHSETPARAMS",ERR,ERROR,*999)

    CALL SNESLineSearchSetParams(SNES_%SNES_,ALPHA,MAXSTEP,STEPTOL,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESLineSearchSetParams",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESLINESEARCHSETPARAMS")
    RETURN
999 CALL ERRORS("PETSC_SNESLINESEARCHSETPARAMS",ERR,ERROR)
    CALL EXITS("PETSC_SNESLINESEARCHSETPARAMS")
    RETURN 1
  END SUBROUTINE PETSC_SNESLINESEARCHSETPARAMS
#endif  
    
  !
  !================================================================================================================================
  !
  
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )
  !>Buffer routine to the PETSc SNESLineSearchSetType routine.
  SUBROUTINE Petsc_SnesLineSearchSetType(lineSearch,lineSearchType,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the line search type for
    SNESLineSearchType, INTENT(IN) :: lineSearchType !<The line search type to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("Petsc_SnesLineSearchSetType",err,error,*999)

    CALL SNESLineSearchSetType(LINESEARCH%snesLineSearch,lineSearchType,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetType",err,error,*999)
    ENDIF

    CALL Exits("Petsc_SnesLineSearchSetType")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchSetType",err,error)
    CALL Exits("Petsc_SnesLineSearchSetType")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchBTSetAlpha routine.
  SUBROUTINE Petsc_SnesLineSearchBTSetAlpha(lineSearch,alpha,err,error,*)

    !Argument variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES back-tracking LineSearch to set alpha for
    REAL(DP), INTENT(IN) :: alpha !<The alpha descent parameter to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("Petsc_SnesLineSearchBTSetAlpha",err,error,*999)

    CALL SNESLineSearchBTSetAlpha(lineSearch%snesLineSearch,alpha,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      END IF
      CALL FlagError("PETSc error in SNESLineSearchBTSetAlpha",err,error,*999)
    END IF

    CALL Exits("Petsc_SnesLineSearchBTSetAlpha")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchBTSetAlpha",err,error)
    CALL Exits("Petsc_SnesLineSearchBTSetAlpha")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchBTSetAlpha

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchSetTolerances routine.
  SUBROUTINE Petsc_SnesLineSearchSetTolerances(lineSearch,steptol,maxstep,rtol,atol,ltol,maxIt,err,error,*)

    !Argument variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set tolerances for
    REAL(DP), INTENT(IN) :: steptol !<The minimum steplength
    REAL(DP), INTENT(IN) :: maxstep !<The maximum steplength
    REAL(DP), INTENT(IN) :: rtol !<The relative tolerance for iterative line searches
    REAL(DP), INTENT(IN) :: atol !<The absolute tolerance for iterative line searches
    REAL(DP), INTENT(IN) :: ltol !<The change in lambda tolerance for iterative line searches
    INTEGER(INTG), INTENT(IN) :: maxIt !<The maximum number of iterations of the line search
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("Petsc_SnesLineSearchSetTolerances",err,error,*999)

    CALL SNESLineSearchSetTolerances(lineSearch%snesLineSearch, &
      & steptol,maxstep,rtol,atol,ltol,maxIt,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      END IF
      CALL FlagError("PETSc error in SNESLineSearchSetTolerances",err,error,*999)
    END IF

    CALL Exits("Petsc_SnesLineSearchSetTolerances")
    RETURN
999 CALL Errors("Petsc_SnesLineSearchSetTolerances",err,error)
    CALL Exits("Petsc_SnesLineSearchSetTolerances")
    RETURN 1

  END SUBROUTINE Petsc_SnesLineSearchSetTolerances
#endif

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESMonitorSet routine.
  SUBROUTINE PETSC_SNESMONITORSET(SNES_,MFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set from the command line options
    EXTERNAL :: MFUNCTION !<The external monitor function to set
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESMONITORSET",ERR,ERROR,*999)

    CALL SNESMonitorSet(SNES_%SNES_,MFUNCTION,CTX,PETSC_NULL_FUNCTION,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESMonitorSet",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESMONITORSET")
    RETURN
999 CALL ERRORS("PETSC_SNESMONITORSET",ERR,ERROR)
    CALL EXITS("PETSC_SNESMONITORSET")
    RETURN 1
  END SUBROUTINE PETSC_SNESMONITORSET
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFromOptions routine.
  SUBROUTINE PETSC_SNESSETFROMOPTIONS(SNES_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set from the command line options
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETFROMOPTIONS",ERR,ERROR,*999)

    CALL SNESSetFromOptions(SNES_%SNES_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetFromOptions",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETFROMOPTIONS")
    RETURN
999 CALL ERRORS("PETSC_SNESSETFROMOPTIONS",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETFROMOPTIONS")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunction routine.
  SUBROUTINE PETSC_SNESSETFUNCTION(SNES_,F,FFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: F !<The residual vector
    EXTERNAL FFUNCTION !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETFUNCTION",ERR,ERROR,*999)

    CALL SNESSetFunction(SNES_%SNES_,F%VEC,FFUNCTION,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetFunction",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETFUNCTION")
    RETURN
999 CALL ERRORS("PETSC_SNESSETFUNCTION",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETFUNCTION")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETFUNCTION

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunction routine.
  SUBROUTINE PETSC_SNESSETCONVERGENCETEST(SNES_,CFUNCTION,CTX,ERR,ERROR,*)
    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
    EXTERNAL CFUNCTION !<The external function to call (OpenCMISS subroutine to calculate convergence
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the convergence test function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETCONVERGENCETEST",ERR,ERROR,*999)

    CALL SNESSetConvergenceTest(SNES_%SNES_,CFUNCTION,CTX,PETSC_NULL_FUNCTION,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetConvergenceTest",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETCONVERGENCETEST")
    RETURN
999 CALL ERRORS("PETSC_SNESSETCONVERGENCETEST",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETCONVERGENCETEST")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETCONVERGENCETEST


  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for MatFDColoring contexts.
  SUBROUTINE PETSC_SNESSETJACOBIAN_MATFDCOLORING(SNES_,A,B,JFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The Jacobian matrix
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
    EXTERNAL JFUNCTION !<The external function to call
    TYPE(PETSC_MATFDCOLORING_TYPE) :: CTX !<The MatFDColoring data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETJACOBIAN_MATFDCOLORING",ERR,ERROR,*999)

    CALL SNESSetJacobianBuffer(SNES_,A,B,JFUNCTION,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetJacobian",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETJACOBIAN_MATFDCOLORING")
    RETURN
999 CALL ERRORS("PETSC_SNESSETJACOBIAN_MATFDCOLORING",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETJACOBIAN_MATFDCOLORING")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETJACOBIAN_MATFDCOLORING
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESGETJACOBIAN_SOLVER(SNES_,A,B,JFUNCTION,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The Jacobian matrix
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
    EXTERNAL JFUNCTION !<The external function to call
!     TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESGETJACOBIAN_SOLVER",ERR,ERROR,*999)

    CALL SNESGetJacobian(SNES_%SNES_,A%MAT,B%MAT,JFUNCTION,PETSC_NULL_INTEGER,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetJacobian",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESGETJACOBIAN_SOLVER")
    RETURN
999 CALL ERRORS("PETSC_SNESGETJACOBIAN_SOLVER",ERR,ERROR)
    CALL EXITS("PETSC_SNESGETJACOBIAN_SOLVER")
    RETURN 1
  END SUBROUTINE PETSC_SNESGETJACOBIAN_SOLVER

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESGETJACOBIAN_SPECIAL(SNES_,A,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The Jacobian matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETJACOBIAN_SPECIAL",ERR,ERROR,*999)

    CALL SNESGetJacobian(SNES_%SNES_,A%MAT,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ERR)

    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetJacobian",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESGETJACOBIAN_SPECIAL")
    RETURN
999 CALL ERRORS("PETSC_SNESGETJACOBIAN_SPECIAL",ERR,ERROR)
    CALL EXITS("PETSC_SNESGETJACOBIAN_SPECIAL")
    RETURN 1
  END SUBROUTINE PETSC_SNESGETJACOBIAN_SPECIAL

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESSETJACOBIAN_SOLVER(SNES_,A,B,JFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The Jacobian matrix
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
    EXTERNAL JFUNCTION !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETJACOBIAN_SOLVER",ERR,ERROR,*999)

    CALL SNESSetJacobian(SNES_%SNES_,A%MAT,B%MAT,JFUNCTION,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetJacobian",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETJACOBIAN_SOLVER")
    RETURN
999 CALL ERRORS("PETSC_SNESSETJACOBIAN_SOLVER",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETJACOBIAN_SOLVER")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETJACOBIAN_SOLVER

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESComputeJacobianColor routine for solver contexts.
  SUBROUTINE PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR(SNES_,X,J,B,FLAG,CTX,ERR,ERROR,*)

    !Argument variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The PETSc SNES
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: J !<The PETSc J Mat
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The PETSc B Mat
    INTEGER(INTG) :: FLAG !<The PETSC MatStructure flag
    TYPE(PETSC_MATFDCOLORING_TYPE), POINTER :: CTX !<The passed through context
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL Enters("PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR",err,error,*999)
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 3 )
    CALL SNESComputeJacobianDefaultColor(SNES_%SNES_,X%VEC,J%MAT,B%MAT,FLAG,CTX,ERR)
#else
    CALL SNESDefaultComputeJacobianColor(SNES_%SNES_,X%VEC,J%MAT,B%MAT,FLAG,CTX,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetJacobian",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR")
    RETURN
999 CALL ERRORS("PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR",ERR,ERROR)
    RETURN
  END SUBROUTINE PETSC_SNESDEFAULTCOMPUTEJACOBIANCOLOR

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESComputeJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESDEFAULTCOMPUTEJACOBIAN(SNES_,X,J,B,FLAG,CTX,ERR,ERROR,*)

    !Argument variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The PETSc SNES
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: J !<The PETSc J Mat
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The PETSc B Mat
    INTEGER(INTG) :: FLAG !<The PETSC MatStructure flag
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("PETSC_SNESDEFAULTCOMPUTEJACOBIAN",ERR,ERROR,*999)
    
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 3 )
    CALL SNESComputeJacobianDefault(SNES_%SNES_,X%VEC,J%MAT,B%MAT,FLAG,CTX,ERR)
#else
    CALL SNESDefaultComputeJacobian(SNES_%SNES_,X%VEC,J%MAT,B%MAT,FLAG,CTX,ERR)
#endif
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetJacobian",ERR,ERROR,*999)
    ENDIF

    RETURN
999 CALL ERRORS("PETSC_SNESDEFAULTCOMPUTEJACOBIAN",ERR,ERROR)
    RETURN
  END SUBROUTINE PETSC_SNESDEFAULTCOMPUTEJACOBIAN

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetKSP routine.
  SUBROUTINE PETSC_SNESSETKSP(SNES_,KSP_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the KSP for
    TYPE(PETSC_KSP_TYPE), INTENT(INOUT) :: KSP_ !<The KSP to be associated with the SNES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETKSP",ERR,ERROR,*999)

    CALL SNESSetKSP(SNES_%SNES_,KSP_%KSP_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetKSP",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETKSP")
    RETURN
999 CALL ERRORS("PETSC_SNESSETKSP",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETKSP")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETKSP

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetTolerances routine.
  SUBROUTINE PETSC_SNESSETTOLERANCES(SNES_,ABSTOL,RTOL,STOL,MAXIT,MAXF,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the tolerances for
    REAL(DP), INTENT(IN) :: ABSTOL !<The absolute convergence tolerance
    REAL(DP), INTENT(IN) :: RTOL !<The relative convergence tolerance
    REAL(DP), INTENT(IN) :: STOL !<The convergence tolerance for the change in the solution between steps
    INTEGER(INTG), INTENT(IN) :: MAXIT !<The maximum number of iterations
    INTEGER(INTG), INTENT(IN) :: MAXF !<The maximum number of function evaluations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETTOLERANCES",ERR,ERROR,*999)

    CALL SNESSetTolerances(SNES_%SNES_,ABSTOL,RTOL,STOL,MAXIT,MAXF,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetTolerances",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETTOLERANCES")
    RETURN
999 CALL ERRORS("PETSC_SNESSETTOLERANCES",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETTOLERANCES")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETTOLERANCES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetTrustRegionTolerance routine.
  SUBROUTINE PETSC_SNESSETTRUSTREGIONTOLERANCE(SNES_,TRTOL,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the tolerances for
    REAL(DP), INTENT(IN) :: TRTOL !<The trust region tolerance
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETTRUSTREGIONTOLERANCE",ERR,ERROR,*999)

    CALL SNESSetTrustRegionTolerance(SNES_%SNES_,TRTOL,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetTrustRegionTolerance",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETTRUSTREGIONTOLERANCE")
    RETURN
999 CALL ERRORS("PETSC_SNESSETTRUSTREGIONTOLERANCE",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETTRUSTREGIONTOLERANCE")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETTRUSTREGIONTOLERANCE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetType routine.
  SUBROUTINE PETSC_SNESSETTYPE(SNES_,METHOD,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the type for
    SNESType, INTENT(IN) :: METHOD !<The SNES type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETTYPE",ERR,ERROR,*999)

    CALL SNESSetType(SNES_%SNES_,METHOD,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETTYPE")
    RETURN
999 CALL ERRORS("PETSC_SNESSETTYPE",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETTYPE")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETTYPE

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchGetVecs routine.
  SUBROUTINE Petsc_SnesLineSearchGetVecs(lineSearch,x,f,y,w,g,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The PetcsSnesLineSearch to get the vectors from the SNESLineSearch
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: x !<The The old solution 
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: f !<The old function 
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: y !<The search direction 
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: w !<The new solution 
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: g !<The new function 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("Petsc_SnesLineSearchGetVecs",err,error,*999)

    CALL SNESLineSearchGetVecs(lineSearch%snesLineSearch,x%VEC,f%VEC,y%VEC,w%VEC,g%VEC,err)
    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESLineSearchGetVecs",err,error,*999)
    ENDIF
    
    CALL EXITS("Petsc_SnesLineSearchGetVecs")
    RETURN
999 CALL ERRORS("Petsc_SnesLineSearchGetVecs",err,error)
    CALL EXITS("Petsc_SnesLineSearchGetVecs")
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchGetVecs

  !
  !================================================================================================================================
  !
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR <= 4 )
  !>Buffer routine to the PETSc SNESSetNormType routine.
  SUBROUTINE PETSC_SNESSETNORMTYPE(SNES_,NORMTYPE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the norm type for
    INTEGER(INTG), INTENT(IN) :: NORMTYPE !<The norm type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSETNORMTYPE",ERR,ERROR,*999)

    CALL SNESSetNormType(SNES_%SNES_,NORMTYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSetNormType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSETNORMTYPE")
    RETURN
999 CALL ERRORS("PETSC_SNESSETNORMTYPE",ERR,ERROR)
    CALL EXITS("PETSC_SNESSETNORMTYPE")
    RETURN 1
  END SUBROUTINE PETSC_SNESSETNORMTYPE
#endif
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSolve routine.
  SUBROUTINE PETSC_SNESSOLVE(SNES_,B,X,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to solve
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: B !<The constant part of the equation
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_SNESSOLVE",ERR,ERROR,*999)

    CALL SNESSolve(SNES_%SNES_,B%VEC,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESSolve",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_SNESSOLVE")
    RETURN
999 CALL ERRORS("PETSC_SNESSOLVE",ERR,ERROR)
    CALL EXITS("PETSC_SNESSOLVE")
    RETURN 1
  END SUBROUTINE PETSC_SNESSOLVE

  !
  !================================================================================================================================
  !

  !Finalise the PETSc TS structure and destroy the TS
  SUBROUTINE PETSC_TSFINALISE(TS_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSFINALISE",ERR,ERROR,*999)

    IF(TS_%TS_/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_TSDESTROY(TS_,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSFINALISE")
    RETURN
999 CALL ERRORS("PETSC_TSFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_TSFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_TSFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc TS structure
  SUBROUTINE PETSC_TSINITIALISE(TS_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSINITIALISE",ERR,ERROR,*999)

    TS_%TS_=PETSC_NULL_INTEGER
     
    CALL EXITS("PETSC_TSINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_TSINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_TSINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_TSINITIALISE
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSCreate routine.
  SUBROUTINE PETSC_TSCREATE(COMMUNICATOR,TS_,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(INOUT) :: COMMUNICATOR !<The MPI communicator for the TS creation
    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<On exit, the TS information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSCREATE",ERR,ERROR,*999)

    CALL TSCreate(COMMUNICATOR,TS_%TS_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSCREATE")
    RETURN
999 CALL ERRORS("PETSC_TSCREATE",ERR,ERROR)
    CALL EXITS("PETSC_TSCREATE")
    RETURN 1
  END SUBROUTINE PETSC_TSCREATE
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSDestroy routine.
  SUBROUTINE PETSC_TSDESTROY(TS_,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSDESTROY",ERR,ERROR,*999)

    CALL TSDestroy(TS_%TS_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSDestroy",ERR,ERROR,*999)
    ENDIF
    TS_%TS_=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_TSDESTROY")
    RETURN
999 CALL ERRORS("PETSC_TSDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_TSDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_TSDESTROY
    
  !
  !================================================================================================================================
  !
    
!  !>Buffer routine to the PETSc TSGetApplicationContext routine.
!  SUBROUTINE PETSC_TSGETAPPLICATIONCONTEXT(TS_,USERP,ERR,ERROR,*)

!    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to get the application context from
!    TYPE(SOLVER_TYPE), POINTER :: USERP !<On exit, a pointer to the user application context
!    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
!    !Local Variables

!    CALL ENTERS("PETSC_TSSETAPPLICATIONCONTEXT",ERR,ERROR,*999)

!    IF(ASSOCIATED(USERP)) THEN
!      CALL FLAG_ERROR("User application pointer is already associated.",ERR,ERROR,*999)
!    ELSE
!      NULLIFY(USERP)
!      CALL TSGetApplicationContext(TS_%TS_,USERP,ERR)
!      IF(ERR/=0) THEN
!        IF(PETSC_HANDLE_ERROR) THEN
!          CHKERRQ(ERR)
!        ENDIF
!        CALL FLAG_ERROR("PETSc error in TSGetApplicationContext",ERR,ERROR,*999)
!      ENDIF
!    ENDIF
    
!    CALL EXITS("PETSC_TSGETAPPLICATIONCONTEXT")
!    RETURN
!999 CALL ERRORS("PETSC_TSGETAPPLICATIONCONTEXT",ERR,ERROR)
!    CALL EXITS("PETSC_TSGETAPPLICATIONCONTEXT")
!    RETURN 1
!  END SUBROUTINE PETSC_TSGETAPPLICATIONCONTEXT
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TSMonitorSet routine.
  SUBROUTINE PETSC_TSMONITORSET(TS_,MFUNCTION,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the monitor for
    EXTERNAL :: MFUNCTION !<The external monitor function to set
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSMONITORSET",ERR,ERROR,*999)

    CALL TSMonitorSet(TS_%TS_,MFUNCTION,CTX,PETSC_NULL_FUNCTION,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSMonitorSet",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSMONITORSET")
    RETURN
999 CALL ERRORS("PETSC_TSMONITORSET",ERR,ERROR)
    CALL EXITS("PETSC_TSMONITORSET")
    RETURN 1
  END SUBROUTINE PETSC_TSMONITORSET
    
  !
  !================================================================================================================================
  !
    
!  !>Buffer routine to the PETSc TSSetApplicationContext routine.
!  SUBROUTINE PETSC_TSSETAPPLICATIONCONTEXT(TS_,USERP,ERR,ERROR,*)

!    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the application context for
!    TYPE(SOLVER_TYPE), POINTER :: USERP !<A pointer to the user application context
!    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
!    !Local Variables

!    CALL ENTERS("PETSC_TSSETAPPLICATIONCONTEXT",ERR,ERROR,*999)

!    CALL TSSetApplicationContext(TS_%TS_,USERP,ERR)
!    IF(ERR/=0) THEN
!      IF(PETSC_HANDLE_ERROR) THEN
!        CHKERRQ(ERR)
!      ENDIF
!      CALL FLAG_ERROR("PETSc error in TSSetApplicationContext",ERR,ERROR,*999)
!    ENDIF
    
!    CALL EXITS("PETSC_TSSETAPPLICATIONCONTEXT")
!    RETURN
!999 CALL ERRORS("PETSC_TSSETAPPLICATIONCONTEXT",ERR,ERROR)
!    CALL EXITS("PETSC_TSSETAPPLICATIONCONTEXT")
!    RETURN 1
!  END SUBROUTINE PETSC_TSSETAPPLICATIONCONTEXT
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetDuration routine.
  SUBROUTINE PETSC_TSSETDURATION(TS_,MAX_STEPS,MAX_TIME,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set from the options
    INTEGER(INTG), INTENT(IN) :: MAX_STEPS !<The maximum number of steps to use
    REAL(DP), INTENT(IN) :: MAX_TIME !<The maximum time to iteration to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETDURATION",ERR,ERROR,*999)

    CALL TSSetDuration(TS_%TS_,MAX_STEPS,MAX_TIME,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetDuration",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETDURATION")
    RETURN
999 CALL ERRORS("PETSC_TSSETDURATION",ERR,ERROR)
    CALL EXITS("PETSC_TSSETDURATION")
    RETURN 1
  END SUBROUTINE PETSC_TSSETDURATION
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetFromOptions routine.
  SUBROUTINE PETSC_TSSETFROMOPTIONS(TS_,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set from the options
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETFROMOPTIONS",ERR,ERROR,*999)

    CALL TSSetFromOptions(TS_%TS_,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetFromOptions",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETFROMOPTIONS")
    RETURN
999 CALL ERRORS("PETSC_TSSETFROMOPTIONS",ERR,ERROR)
    CALL EXITS("PETSC_TSSETFROMOPTIONS")
    RETURN 1
  END SUBROUTINE PETSC_TSSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetInitialTimeStep routine.
  SUBROUTINE PETSC_TSSETINITIALTIMESTEP(TS_,INITIAL_TIME,TIME_STEP,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the initial time step for
    REAL(DP), INTENT(IN) :: INITIAL_TIME !<The initial time to set
    REAL(DP), INTENT(IN) :: TIME_STEP !<The time step to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETINITIALTIMESTEP",ERR,ERROR,*999)

    CALL TSSetInitialTimeStep(TS_%TS_,INITIAL_TIME,TIME_STEP,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetInitialTimeStep",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETINITIALTIMESTEP")
    RETURN
999 CALL ERRORS("PETSC_TSSETINITIALTIMESTEP",ERR,ERROR)
    CALL EXITS("PETSC_TSSETINITIALTIMESTEP")
    RETURN 1
  END SUBROUTINE PETSC_TSSETINITIALTIMESTEP

  !
  !================================================================================================================================
  !
    
#if ( PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2 )
  !>Buffer routine to the PETSc TSSetMatrices routine.
  SUBROUTINE PETSC_TSSETMATRICES(TS_,ARHS,RHSFUNCTION,ALHS,LHSFUNCTION,FLAG,CTX,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the problem type for
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: ARHS !<The RHS matrix
    EXTERNAL RHSFUNCTION !<The external RHS matrix evaluation function to call
    TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: ALHS !<The LHS matrix
    EXTERNAL LHSFUNCTION !<The external LHS matrix evaluation function to call
    INTEGER(INTG), INTENT(IN) :: FLAG !<The matrices structure flag
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the matrix evaluations functions
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETMATRICES",ERR,ERROR,*999)

    CALL TSSetMatrices(TS_%TS_,ARHS%MAT,RHSFUNCTION,ALHS%MAT,LHSFUNCTION,FLAG,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetMatrices",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETMATRICES")
    RETURN
999 CALL ERRORS("PETSC_TSSETMATRICES",ERR,ERROR)
    CALL EXITS("PETSC_TSSETMATRICES")
    RETURN 1
  END SUBROUTINE PETSC_TSSETMATRICES
    
  !
  !================================================================================================================================
  !
#endif
  
  !>Buffer routine to the PETSc TSSetProblemType routine.
  SUBROUTINE PETSC_TSSETPROBLEMTYPE(TS_,PROB_TYPE,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the problem type for
    INTEGER(INTG), INTENT(IN) :: PROB_TYPE !<The problem type to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETPROBLEMTYPE",ERR,ERROR,*999)

    CALL TSSetProblemType(TS_%TS_,PROB_TYPE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetProblemType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETPROBLEMTYPE")
    RETURN
999 CALL ERRORS("PETSC_TSSETPROBLEMTYPE",ERR,ERROR)
    CALL EXITS("PETSC_TSSETPROBLEMTYPE")
    RETURN 1
  END SUBROUTINE PETSC_TSSETPROBLEMTYPE
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetRHSFunction routine.
  SUBROUTINE PETSC_TSSETRHSFUNCTION(TS_,RHSFUNCTION,CTX,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the problem type for
    EXTERNAL RHSFUNCTION !<The external RHS function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETRHSFUNCTION",ERR,ERROR,*999)

    CALL TSSetRHSFunction(TS_%TS_,RHSFUNCTION,CTX,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetRHSFunction",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETRHSFUNCTION")
    RETURN
999 CALL ERRORS("PETSC_TSSETRHSFUNCTION",ERR,ERROR)
    CALL EXITS("PETSC_TSSETRHSFUNCTION")
    RETURN 1
  END SUBROUTINE PETSC_TSSETRHSFUNCTION
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetTimeStep routine.
  SUBROUTINE PETSC_TSSETTIMESTEP(TS_,TIME_STEP,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the time step for
    REAL(DP), INTENT(IN) :: TIME_STEP !<The time step to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETTIMESTEP",ERR,ERROR,*999)

    CALL TSSetTimeStep(TS_%TS_,TIME_STEP,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetTimeStep",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETTIMESTEP")
    RETURN
999 CALL ERRORS("PETSC_TSSETTIMESTEP",ERR,ERROR)
    CALL EXITS("PETSC_TSSETTIMESTEP")
    RETURN 1
  END SUBROUTINE PETSC_TSSETTIMESTEP
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetType routine.
  SUBROUTINE PETSC_TSSETTYPE(TS_,METHOD,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to set the type for
    TSType, INTENT(IN) :: METHOD !<The time stepping method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSETTYPE",ERR,ERROR,*999)

    CALL TSSetType(TS_%TS_,METHOD,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSetType",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSETTYPE")
    RETURN
999 CALL ERRORS("PETSC_TSSETTYPE",ERR,ERROR)
    CALL EXITS("PETSC_TSSETTYPE")
    RETURN 1
  END SUBROUTINE PETSC_TSSETTYPE
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSolve routine.
  SUBROUTINE PETSC_TSSOLVE(TS_,X,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to solve
    TYPE(PETSC_VEC_TYPE), INTENT(IN) :: X !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSOLVE",ERR,ERROR,*999)

    CALL TSSolve(TS_%TS_,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSSolve",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSOLVE")
    RETURN
999 CALL ERRORS("PETSC_TSSOLVE",ERR,ERROR)
    CALL EXITS("PETSC_TSSOLVE")
    RETURN 1
  END SUBROUTINE PETSC_TSSOLVE
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSStep routine.
  SUBROUTINE PETSC_TSSTEP(TS_,STEPS,PTIME,ERR,ERROR,*)

    TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS_ !<The TS to step
    INTEGER(INTG), INTENT(IN) :: STEPS !<The number of iterations until termination
    REAL(DP), INTENT(IN) :: PTIME !<The time until termination
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_TSSTEP",ERR,ERROR,*999)

    CALL TSStep(TS_%TS_,STEPS,PTIME,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in TSStep",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_TSSTEP")
    RETURN
999 CALL ERRORS("PETSC_TSSTEP",ERR,ERROR)
    CALL EXITS("PETSC_TSSTEP")
    RETURN 1
  END SUBROUTINE PETSC_TSSTEP
  
  !
  !================================================================================================================================
  !

  !Finalise the PETSc Vec structure and destroy the KSP
  SUBROUTINE PETSC_VECFINALISE(VEC_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: VEC_ !<The Vec to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECFINALISE",ERR,ERROR,*999)

    IF(VEC_%VEC/=PETSC_NULL_INTEGER) THEN
      CALL PETSC_VECDESTROY(VEC_,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECFINALISE")
    RETURN
999 CALL ERRORS("PETSC_VECFINALISE",ERR,ERROR)
    CALL EXITS("PETSC_VECFINALISE")
    RETURN 1
  END SUBROUTINE PETSC_VECFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Vec structure
  SUBROUTINE PETSC_VECINITIALISE(VEC_,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: VEC_ !<The Vec to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECINITIALISE",ERR,ERROR,*999)

    VEC_%VEC=PETSC_NULL_INTEGER
    !VEC_%VEC_DATA(1)=0
    !VEC_%VEC_OFFSET=0
    
    CALL EXITS("PETSC_VECINITIALISE")
    RETURN
999 CALL ERRORS("PETSC_VECINITIALISE",ERR,ERROR)
    CALL EXITS("PETSC_VECINITIALISE")
    RETURN 1
  END SUBROUTINE PETSC_VECINITIALISE
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecAssemblyBegin routine.
  SUBROUTINE PETSC_VECASSEMBLYBEGIN(X,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to begin the assembly of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECASSEMBLYBEGIN",ERR,ERROR,*999)

    CALL VecAssemblyBegin(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecAssemblyBegin",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECASSEMBLYBEGIN")
    RETURN
999 CALL ERRORS("PETSC_VECASSEMBLYBEGIN",ERR,ERROR)
    CALL EXITS("PETSC_VECASSEMBLYBEGIN")
    RETURN 1
  END SUBROUTINE PETSC_VECASSEMBLYBEGIN
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecAssemblyEnd routine.
  SUBROUTINE PETSC_VECASSEMBLYEND(X,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to end the assembly of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECASSEMBLYEND",ERR,ERROR,*999)

    CALL VecAssemblyEnd(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecAssemblyEnd",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECASSEMBLYEND")
    RETURN
999 CALL ERRORS("PETSC_VECASSEMBLYEND",ERR,ERROR)
    CALL EXITS("PETSC_VECASSEMBLYEND")
    RETURN 1
  END SUBROUTINE PETSC_VECASSEMBLYEND
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCopy routine.
  SUBROUTINE PETSC_VECCOPY(X,Y,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to copy from
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: Y !<The vector to copy to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCOPY",ERR,ERROR,*999)

    CALL VecCopy(X%VEC,Y%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCopy",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCOPY")
    RETURN
999 CALL ERRORS("PETSC_VECCOPY",ERR,ERROR)
    CALL EXITS("PETSC_VECCOPY")
    RETURN 1
    
  END SUBROUTINE PETSC_VECCOPY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreate routine.
  SUBROUTINE PETSC_VECCREATE(COMMUNICATOR,X,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATE",ERR,ERROR,*999)

    CALL VecCreate(COMMUNICATOR,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATE")
    RETURN
999 CALL ERRORS("PETSC_VECCREATE",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATE")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateGhost routine.
  SUBROUTINE PETSC_VECCREATEGHOST(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,X,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    INTEGER(INTG), INTENT(IN) :: NUMBER_GHOST !<The number of ghost elements
    INTEGER(INTG), INTENT(IN) :: GHOSTS(*) !<The global location of the each ghost element
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEGHOST",ERR,ERROR,*999)

    CALL VecCreateGhost(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateGhost",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEGHOST")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEGHOST",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEGHOST")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEGHOST
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateGhostWithArray routine.
  SUBROUTINE PETSC_VECCREATEGHOSTWITHARRAY(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,ARRAY,X,ERR,ERROR,*)

   !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    INTEGER(INTG), INTENT(IN) :: NUMBER_GHOST !<The number of ghost elements
    INTEGER(INTG), INTENT(IN) :: GHOSTS(*) !<The global location of the each ghost element
    REAL(DP), INTENT(OUT) :: ARRAY(*) !<The preallocated array of matrix data
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEGHOSTWITHARRAY",ERR,ERROR,*999)

    CALL VecCreateGhostWithArray(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,ARRAY,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateGhostWithArray",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEGHOSTWITHARRAY")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEGHOSTWITHARRAY",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEGHOSTWITHARRAY")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEGHOSTWITHARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateMPI routine.
  SUBROUTINE PETSC_VECCREATEMPI(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,X,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEMPI",ERR,ERROR,*999)

    CALL VecCreateMPI(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateMPI",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEMPI")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEMPI",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEMPI")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEMPI
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateMPIWithArray routine.
  SUBROUTINE PETSC_VECCREATEMPIWITHARRAY(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,ARRAY,X,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    REAL(DP), INTENT(OUT) :: ARRAY(*) !<The preallocated array for the vector data
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATEMPIWITHARRAY",ERR,ERROR,*999)

    CALL VecCreateMPIWithArray(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,ARRAY,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateMPIWithArray",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATEMPIWITHARRAY")
    RETURN
999 CALL ERRORS("PETSC_VECCREATEMPIWITHARRAY",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATEMPIWITHARRAY")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEMPIWITHARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateSeq routine.
  SUBROUTINE PETSC_VECCREATESEQ(COMMUNICATOR,SIZE,X,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: SIZE !<The size of the vector
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATESEQ",ERR,ERROR,*999)

    CALL VecCreateSeq(COMMUNICATOR,SIZE,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateSeq",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATESEQ")
    RETURN
999 CALL ERRORS("PETSC_VECCREATESEQ",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATESEQ")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATESEQ
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateSeqWithArray routine.
  SUBROUTINE PETSC_VECCREATESEQWITHARRAY(COMMUNICATOR,SIZE,ARRAY,X,ERR,ERROR,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: SIZE !<The size of the vector
    REAL(DP), INTENT(OUT) :: ARRAY(*) !<The preallocated array for the vector data
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECCREATESEQWITHARRAY",ERR,ERROR,*999)

    CALL VecCreateSeqWithArray(COMMUNICATOR,SIZE,ARRAY,X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecCreateSeqWithArray",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECCREATESEQWITHARRAY")
    RETURN
999 CALL ERRORS("PETSC_VECCREATESEQWITHARRAY",ERR,ERROR)
    CALL EXITS("PETSC_VECCREATESEQWITHARRAY")
    RETURN 1
  END SUBROUTINE PETSC_VECCREATESEQWITHARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDestroy routine.
  SUBROUTINE PETSC_VECDESTROY(X,ERR,ERROR,*)

   !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECDESTROY",ERR,ERROR,*999)

    CALL VecDestroy(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecDestroy",ERR,ERROR,*999)
    ENDIF
    X%VEC=PETSC_NULL_INTEGER
    
    CALL EXITS("PETSC_VECDESTROY")
    RETURN
999 CALL ERRORS("PETSC_VECDESTROY",ERR,ERROR)
    CALL EXITS("PETSC_VECDESTROY")
    RETURN 1
  END SUBROUTINE PETSC_VECDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDuplicate routine.
  SUBROUTINE PETSC_VECDUPLICATE(OLD,NEW,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: OLD !<The vector to duplicate
    TYPE(PETSC_VEC_TYPE), INTENT(OUT) :: NEW !<On exit, the new duplicated vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECDUPLICATE",ERR,ERROR,*999)

    CALL VecDuplicate(OLD%VEC,NEW%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecDuplicate",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECDUPLICATE")
    RETURN
999 CALL ERRORS("PETSC_VECDUPLICATE",ERR,ERROR)
    CALL EXITS("PETSC_VECDUPLICATE")
    RETURN 1
  END SUBROUTINE PETSC_VECDUPLICATE
    
  !
  !================================================================================================================================
  !

!  !>Buffer routine to the PETSc VecGetArray routine.
!  SUBROUTINE PETSC_VECGETARRAY(X,ARRAY,ERR,ERROR,*)

!    !Argument Variables
!    TYPE(PETSC_VEC_TYPE), INTENT(INOUT), TARGET :: X !<The vector to get the array of
!    REAL(DP), POINTER :: ARRAY(:) !<On exit, a pointer to the array of the vector
!    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
!    !Local Variables

!    CALL ENTERS("PETSC_VECGETARRAY",ERR,ERROR,*999)

!    IF(ASSOCIATED(ARRAY)) THEN
!      CALL FLAG_ERROR("Array is already associated",ERR,ERROR,*999)
!    ELSE
!      CALL VecGetArray(X%VEC,X%VEC_DATA,X%VEC_OFFSET,ERR)
!      IF(ERR/=0) THEN
!        IF(PETSC_HANDLE_ERROR) THEN
!          CHKERRQ(ERR)
!        ENDIF
!        CALL FLAG_ERROR("PETSc error in VecGetArray",ERR,ERROR,*999)
!      ENDIF
!      ARRAY=>X%VEC_DATA(X%VEC_OFFSET:)
!    ENDIF
    
!    CALL EXITS("PETSC_VECGETARRAY")
!    RETURN
!999 CALL ERRORS("PETSC_VECGETARRAY",ERR,ERROR)
!    CALL EXITS("PETSC_VECGETARRAY")
!    RETURN 1
!  END SUBROUTINE PETSC_VECGETARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetArrayF90 routine.
  SUBROUTINE PETSC_VECGETARRAYF90(X,ARRAY,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT), TARGET :: X !<The vector to get the array of
    REAL(DP), POINTER :: ARRAY(:) !<On exit, a pointer to the array of the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGETARRAYF90",ERR,ERROR,*999)

    IF(ASSOCIATED(ARRAY)) THEN
      CALL FLAG_ERROR("Array is already associated",ERR,ERROR,*999)
    ELSE
      CALL VecGetArrayF90(X%VEC,ARRAY,ERR)
      IF(ERR/=0) THEN
        IF(PETSC_HANDLE_ERROR) THEN
          CHKERRQ(ERR)
        ENDIF
        CALL FLAG_ERROR("PETSc error in VecGetArrayF90",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    CALL EXITS("PETSC_VECGETARRAYF90")
    RETURN
999 CALL ERRORS("PETSC_VECGETARRAYF90",ERR,ERROR)
    CALL EXITS("PETSC_VECGETARRAYF90")
    RETURN 1
  END SUBROUTINE PETSC_VECGETARRAYF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetLocalSize routine.
  SUBROUTINE PETSC_VECGETLOCALSIZE(X,SIZE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to get the local size of
    INTEGER(INTG), INTENT(OUT) :: SIZE !<On exit, the local size of the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGETLOCALSIZE",ERR,ERROR,*999)

    CALL VecGetLocalSize(X%VEC,SIZE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetLocalSize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETLOCALSIZE")
    RETURN
999 CALL ERRORS("PETSC_VECGETLOCALSIZE",ERR,ERROR)
    CALL EXITS("PETSC_VECGETLOCALSIZE")
    RETURN 1
  END SUBROUTINE PETSC_VECGETLOCALSIZE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetOwnershipRange routine.
  SUBROUTINE PETSC_VECGETOWNERSHIPRANGE(X,LOW,HIGH,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to get the ownership range of 
    INTEGER(INTG), INTENT(OUT) :: LOW !<On exit, the low end of the range
    INTEGER(INTG), INTENT(OUT) :: HIGH !<On exit, the high end of the range
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGETOWNERSHIPRANGE",ERR,ERROR,*999)

    CALL VecGetOwnershipRange(X%VEC,LOW,HIGH,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
        ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetOwnershipRange",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETOWNERSHIPRANGE")
    RETURN
999 CALL ERRORS("PETSC_VECGETOWNERSHIPRANGE",ERR,ERROR)
    CALL EXITS("PETSC_VECGETOWNERSHIPRANGE")
    RETURN 1
  END SUBROUTINE PETSC_VECGETOWNERSHIPRANGE

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDot routine.
  SUBROUTINE Petsc_VecDot(x,y,dotProduct,err,error,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(IN) :: x !<The vector x
    TYPE(PETSC_VEC_TYPE), INTENT(IN) :: y !<The vector y
    REAL(DP), INTENT(OUT) :: dotProduct !<The dot product 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL ENTERS("Petsc_VecDot",err,error,*999)
    
    CALL VecDot(x%VEC,y%VEC,dotProduct,err)

    IF(err/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(err)
      ENDIF
      CALL FLAG_ERROR("PETSc error in SNESGetSolutionUpdate",err,error,*999)
    ENDIF

    CALL EXITS("Petsc_VecDot")
    RETURN
999 CALL ERRORS("Petsc_VecDot",err,error)
    CALL EXITS("Petsc_VecDot")
    RETURN 1
  END SUBROUTINE Petsc_VecDot

    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetSize routine.
  SUBROUTINE PETSC_VECGETSIZE(X,SIZE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to get the size of
    INTEGER(INTG), INTENT(OUT) :: SIZE !<On exit, the size of the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGETSIZE",ERR,ERROR,*999)

    CALL VecGetSize(X%VEC,SIZE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetSize",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETSIZE")
    RETURN
999 CALL ERRORS("PETSC_VECGETSIZE",ERR,ERROR)
    CALL EXITS("PETSC_VECGETSIZE")
    RETURN 1
  END SUBROUTINE PETSC_VECGETSIZE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetValues routine.
  SUBROUTINE PETSC_VECGETVALUES(X,N,INDICES,VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the values for
    INTEGER(INTG), INTENT(IN) :: N !<The number of indicies to get
    INTEGER(INTG), INTENT(IN) :: INDICES(*) !<The indices to get
    REAL(DP), INTENT(OUT) :: VALUES(*) !<On return, the values at the specified indicies
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGETVALUES",ERR,ERROR,*999)

    CALL VecGetValues(X%VEC,N,INDICES,VALUES,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGetValues",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGETVALUES")
    RETURN
999 CALL ERRORS("PETSC_VECGETVALUES",ERR,ERROR)
    CALL EXITS("PETSC_VECGETVALUES")
    RETURN 1
  END SUBROUTINE PETSC_VECGETVALUES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostGetLocalForm routine.
  SUBROUTINE PETSC_VECGHOSTGETLOCALFORM(G,L,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: G !<The global form of the vector
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: L !<On exit, the local form of the vector with ghosts
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTGETLOCALFORM",ERR,ERROR,*999)

    CALL VecGhostGetLocalForm(G%VEC,L%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostGetLocalForm",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTGETLOCALFORM")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTGETLOCALFORM",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTGETLOCALFORM")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTGETLOCALFORM
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostRestoreLocalForm routine.
  SUBROUTINE PETSC_VECGHOSTRESTORELOCALFORM(G,L,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: G !<The global form of the vector
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: L !<The local form of the vector
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTRESTORELOCALFORM",ERR,ERROR,*999)

    CALL VecGhostRestoreLocalForm(G%VEC,L%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostRestoreLocalForm",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTRESTORELOCALFORM")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTRESTORELOCALFORM",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTRESTORELOCALFORM")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTRESTORELOCALFORM
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostUpdateBegin routine.
  SUBROUTINE PETSC_VECGHOSTUPDATEBEGIN(X,INSERT_MODE,SCATTER_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to begin the ghost update for
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    ScatterMode, INTENT(IN) :: SCATTER_MODE !<The scatter mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTUPDATEBEGIN",ERR,ERROR,*999)

    CALL VecGhostUpdateBegin(X%VEC,INSERT_MODE,SCATTER_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostUpdateBegin",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTUPDATEBEGIN")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTUPDATEBEGIN",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTUPDATEBEGIN")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTUPDATEBEGIN
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostUpdateEnd routine.
  SUBROUTINE PETSC_VECGHOSTUPDATEEND(X,INSERT_MODE,SCATTER_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to end the ghost update for
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    ScatterMode, INTENT(IN) :: SCATTER_MODE !<The scatter mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECGHOSTUPDATEEND",ERR,ERROR,*999)

    CALL VecGhostUpdateEnd(X%VEC,INSERT_MODE,SCATTER_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecGhostUpdateEnd",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECGHOSTUPDATEEND")
    RETURN
999 CALL ERRORS("PETSC_VECGHOSTUPDATEEND",ERR,ERROR)
    CALL EXITS("PETSC_VECGHOSTUPDATEEND")
    RETURN 1
  END SUBROUTINE PETSC_VECGHOSTUPDATEEND
    
  !
  !================================================================================================================================
  !

!  !>Buffer routine to the PETSc VecRestoreArray routine.
!  SUBROUTINE PETSC_VECRESTOREARRAY(X,ERR,ERROR,*)

!    !Argument Variables
!    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to restore the array of
!    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
!    !Local Variables

!    CALL ENTERS("PETSC_VECRESTOREARRAY",ERR,ERROR,*999)

!    CALL VecRestoreArray(X%VEC,X%VEC_DATA,X%VEC_OFFSET,ERR)
!    IF(ERR/=0) THEN
!      IF(PETSC_HANDLE_ERROR) THEN
!        CHKERRQ(ERR)
!      ENDIF
!      CALL FLAG_ERROR("PETSc error in VecRestoreArray",ERR,ERROR,*999)
!    ENDIF
    
!    CALL EXITS("PETSC_VECRESTOREARRAY")
!    RETURN
!999 CALL ERRORS("PETSC_VECRESTOREARRAY",ERR,ERROR)
!    CALL EXITS("PETSC_VECRESTOREARRAY")
!    RETURN 1
!  END SUBROUTINE PETSC_VECRESTOREARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecRestoreArrayF90 routine.
  SUBROUTINE PETSC_VECRESTOREARRAYF90(X,ARRAY,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to restore the array of
    REAL(DP), POINTER :: ARRAY(:) !<A pointer to the data to restore
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECRESTOREARRAYF90",ERR,ERROR,*999)

    CALL VecRestoreArrayF90(X%VEC,ARRAY,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecRestoreArrayF90",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECRESTOREARRAYF90")
    RETURN
999 CALL ERRORS("PETSC_VECRESTOREARRAYF90",ERR,ERROR)
    CALL EXITS("PETSC_VECRESTOREARRAYF90")
    RETURN 1
  END SUBROUTINE PETSC_VECRESTOREARRAYF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecScale routine.
  SUBROUTINE PETSC_VECSCALE(X,ALPHA,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to scale
    REAL(DP), INTENT(IN) :: ALPHA !<The scaling factor
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSCALE",ERR,ERROR,*999)

    CALL VecScale(X%VEC,ALPHA,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecScale",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSCALE")
    RETURN
999 CALL ERRORS("PETSC_VECSCALE",ERR,ERROR)
    CALL EXITS("PETSC_VECSCALE")
    RETURN 1
    
  END SUBROUTINE PETSC_VECSCALE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSet routine.
  SUBROUTINE PETSC_VECSET(X,VALUE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the value of
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSET",ERR,ERROR,*999)

    CALL VecSet(X%VEC,VALUE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSet",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSET")
    RETURN
999 CALL ERRORS("PETSC_VECSET",ERR,ERROR)
    CALL EXITS("PETSC_VECSET")
    RETURN 1
  END SUBROUTINE PETSC_VECSET
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetFromOptions routine.
  SUBROUTINE PETSC_VECSETFROMOPTIONS(X,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the options for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSETFROMOPTIONS",ERR,ERROR,*999)

    CALL VecSetFromOptions(X%VEC,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetFromOptions",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETFROMOPTIONS")
    RETURN
999 CALL ERRORS("PETSC_VECSETFROMOPTIONS",ERR,ERROR)
    CALL EXITS("PETSC_VECSETFROMOPTIONS")
    RETURN 1
  END SUBROUTINE PETSC_VECSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetLocalToGlobalMapping routine.
  SUBROUTINE PETSC_VECSETLOCALTOGLOBALMAPPING(X,CTX,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the local to global mapping for
    TYPE(PETSC_ISLOCALTOGLOBALMAPPING_TYPE), INTENT(IN) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSETLOCALTOGLOBALMAPPING",ERR,ERROR,*999)

    CALL VecSetLocalToGlobalMapping(X%VEC,CTX%ISLOCALTOGLOBALMAPPING,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetLocalToGlobalMapping",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETLOCALTOGLOBALMAPPING")
    RETURN
999 CALL ERRORS("PETSC_VECSETLOCALTOGLOBALMAPPING",ERR,ERROR)
    CALL EXITS("PETSC_VECSETLOCALTOGLOBALMAPPING")
    RETURN 1
  END SUBROUTINE PETSC_VECSETLOCALTOGLOBALMAPPING
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetValues routine.
  SUBROUTINE PETSC_VECSETVALUES(X,N,INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the values for
    INTEGER(INTG), INTENT(IN) :: N !<The number of indicies
    INTEGER(INTG), INTENT(IN) :: INDICES(*) !<The indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSETVALUES",ERR,ERROR,*999)

    CALL VecSetValues(X%VEC,N,INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetValues",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETVALUES")
    RETURN
999 CALL ERRORS("PETSC_VECSETVALUES",ERR,ERROR)
    CALL EXITS("PETSC_VECSETVALUES")
    RETURN 1
  END SUBROUTINE PETSC_VECSETVALUES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetValuesLocal routine.
  SUBROUTINE PETSC_VECSETVALUESLOCAL(X,N,INDICES,VALUES,INSERT_MODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the values of
    INTEGER(INTG), INTENT(IN) :: N !<The number of indices
    INTEGER(INTG), INTENT(IN) :: INDICES(*) !<The local indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSETVALUESLOCAL",ERR,ERROR,*999)

    CALL VecSetValuesLocal(X%VEC,N,INDICES,VALUES,INSERT_MODE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetValuesLocal",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETVALUESLOCAL")
    RETURN
999 CALL ERRORS("PETSC_VECSETVALUESLOCAL",ERR,ERROR)
    CALL EXITS("PETSC_VECSETVALUESLOCAL")
    RETURN 1
  END SUBROUTINE PETSC_VECSETVALUESLOCAL
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetSizes routine.
  SUBROUTINE PETSC_VECSETSIZES(X,LOCAL_SIZE,GLOBAL_SIZE,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to set the sizes of
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECSETSIZES",ERR,ERROR,*999)

    CALL VecSetSizes(X%VEC,LOCAL_SIZE,GLOBAL_SIZE,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecSetSizes",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECSETSIZES")
    RETURN
999 CALL ERRORS("PETSC_VECSETSIZES",ERR,ERROR)
    CALL EXITS("PETSC_VECSETSIZES")
    RETURN 1
  END SUBROUTINE PETSC_VECSETSIZES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecView routine.
  SUBROUTINE PETSC_VECVIEW(X,V,ERR,ERROR,*)

    !Argument Variables
    TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The vector to view
    PetscViewer, INTENT(IN) :: V !<The viewer
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PETSC_VECVIEW",ERR,ERROR,*999)

    CALL VecView(X%VEC,V,ERR)
    IF(ERR/=0) THEN
      IF(PETSC_HANDLE_ERROR) THEN
        CHKERRQ(ERR)
      ENDIF
      CALL FLAG_ERROR("PETSc error in VecView",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PETSC_VECVIEW")
    RETURN
999 CALL ERRORS("PETSC_VECVIEW",ERR,ERROR)
    CALL EXITS("PETSC_VECVIEW")
    RETURN 1
  END SUBROUTINE PETSC_VECVIEW
    
  !
  !================================================================================================================================
  !

END MODULE CMISS_PETSC
    
!>Buffer routine to the PETSc SNESSetJacobian routine for MatFDColoring contexts. The buffer is required because we want to
!>provide an interface so that we can pass a pointer to the solver for analytic Jacobian's. However, if we provided an interface
!>the Fortran's strong typing rules would not let us pass the matfdcoloring.
SUBROUTINE SNESSetJacobianBuffer(SNES_,A,B,JFUNCTION,CTX,ERR)

  USE CMISS_PETSC_TYPES
  USE KINDS

  IMPLICIT NONE
  
  !Argument Variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES_ !<The SNES to set the function for
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The Jacobian matrix
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
  EXTERNAL JFUNCTION !<The external function to call
  TYPE(PETSC_MATFDCOLORING_TYPE) :: CTX !<The MatFDColoring data to pass to the function
  INTEGER(INTG), INTENT(OUT) :: ERR !<The error code

  CALL SNESSetJacobian(SNES_%SNES_,A%MAT,B%MAT,JFUNCTION,CTX%MATFDCOLORING,ERR)

END SUBROUTINE SNESSetJacobianBuffer
