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
  USE CmissPetscTypes
  USE KINDS
  USE ISO_VARYING_STRING
  USE TYPES
  
#include "macros.h"

  IMPLICIT NONE
 
  PRIVATE

#include "petscversion.h"

#if ( PETSC_VERSION_GE(3,6,0) )
#include "petsc/finclude/petsc.h"
#else
#include "finclude/petsc.h"
#if ( PETSC_VERSION_LT(3,1,0) )
#include "finclude/petscis.h"
#include "finclude/petscksp.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscsnes.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#endif
#if ( PETSC_VERSION_LT(3,3,0) )
#include "finclude/petscts.h"
#endif
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
#if ( PETSC_VERSION_LT(3,4,0) )
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NAN = KSP_DIVERGED_NAN
#else
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NAN = KSP_DIVERGED_NANORINF
#endif
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_INDEFINITE_MAT = KSP_DIVERGED_INDEFINITE_MAT
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_ITERATING = KSP_CONVERGED_ITERATING  
  
  !KSP types
  KSPType, PARAMETER :: PETSC_KSPRICHARDSON = KSPRICHARDSON
#if ( PETSC_VERSION_GE(3,4,0) )
  KSPType, PARAMETER :: PETSC_KSPCHEBYSHEV = KSPCHEBYSHEV
#else
  KSPType, PARAMETER :: PETSC_KSPCHEBYCHEV = KSPCHEBYCHEV
#endif  
  KSPType, PARAMETER :: PETSC_KSPCG = KSPCG
  KSPType, PARAMETER :: PETSC_KSPCGNE = KSPCGNE
#if ( PETSC_VERSION_GE(3,5,0) )
  KSPType, PARAMETER :: PETSC_KSPNASH = KSPNASH
#endif
  KSPType, PARAMETER :: PETSC_KSPSTCG = KSPSTCG
#if ( PETSC_VERSION_GE(3,5,0) )
  KSPType, PARAMETER :: PETSC_KSPGLTR = KSPGLTR
#endif
  KSPType, PARAMETER :: PETSC_KSPGMRES = KSPGMRES
  KSPType, PARAMETER :: PETSC_KSPFGMRES = KSPFGMRES
  KSPType, PARAMETER :: PETSC_KSPLGMRES = KSPLGMRES
#if ( PETSC_VERSION_GE(3,5,0) )
  KSPType, PARAMETER :: PETSC_KSPDGMRES = KSPDGMRES
  KSPType, PARAMETER :: PETSC_KSPPGMRES = KSPPGMRES
#endif
  KSPType, PARAMETER :: PETSC_KSPTCQMR = KSPTCQMR
  KSPType, PARAMETER :: PETSC_KSPBCGS = KSPBCGS
#if ( PETSC_VERSION_GE(3,5,0) )
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
#if ( PETSC_VERSION_GE(3,5,0) )
  KSPType, PARAMETER :: PETSC_KSPPYTHON = KSPPYTHON
  KSPType, PARAMETER :: PETSC_KSPGCR = KSPGCR
#if ( PETSC_VERSION_LE(3,5,4) )
  KSPType, PARAMETER :: PETSC_KSPSPECEST = KSPSPECEST
#endif
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
#if ( PETSC_VERSION_GE(3,5,0) )
  MatOption, PARAMETER :: PETSC_MAT_DUMMY = MAT_DUMMY
#endif  
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_LOWER_TRIANGULAR = MAT_IGNORE_LOWER_TRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_ERROR_LOWER_TRIANGULAR = MAT_ERROR_LOWER_TRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_GETROW_UPPERTRIANGULAR = MAT_GETROW_UPPERTRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR = MAT_UNUSED_NONZERO_LOCATION_ERR
#if ( PETSC_VERSION_GE(3,4,0) )
  MatOption, PARAMETER :: PETSC_NUM_MAT_OPTIONS = MAT_OPTION_MAX
#else
  MatOption, PARAMETER :: PETSC_NUM_MAT_OPTIONS = NUM_MAT_OPTIONS
#endif
#if ( PETSC_VERSION_GE(3,2,0) )
#if ( PETSC_VERSION_LT(3,5,0) )
  MatOption, PARAMETER :: PETSC_MAT_CHECK_COMPRESSED_ROW = MAT_CHECK_COMPRESSED_ROW
#endif  
#else
  MatOption, PARAMETER :: PETSC_MAT_USE_COMPRESSEDROW = MAT_USE_COMPRESSEDROW
#endif  
#if ( PETSC_VERSION_GE(3,2,0) )
  MatOption, PARAMETER :: PETSC_MAT_SPD = MAT_SPD
  MatOption, PARAMETER :: PETSC_MAT_NO_OFF_PROC_ENTRIES = MAT_NO_OFF_PROC_ENTRIES
  MatOption, PARAMETER :: PETSC_MAT_NO_OFF_PROC_ZERO_ROWS = MAT_NO_OFF_PROC_ZERO_ROWS
#endif
  
  !MatStructure types
#if ( PETSC_VERSION_LT(3,5,0) )
  MatStructure, PARAMETER :: PETSC_SAME_PRECONDITIONER = SAME_PRECONDITIONER
#endif  
  MatStructure, PARAMETER :: PETSC_SAME_NONZERO_PATTERN = SAME_NONZERO_PATTERN
  MatStructure, PARAMETER :: PETSC_DIFFERENT_NONZERO_PATTERN = DIFFERENT_NONZERO_PATTERN
  MatStructure, PARAMETER :: PETSC_SUBSET_NONZERO_PATTERN = SUBSET_NONZERO_PATTERN

  !MatColoring types
#if ( PETSC_VERSION_GE(3,2,0) )
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
#if ( PETSC_VERSION_GE(3,5,0) )
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_GREEDY = MATCOLORINGGREEDY
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_JP = MATCOLORINGJP
#endif

  !Norm types
  NormType, PARAMETER :: PETSC_NORM_1 = NORM_1
  NormType, PARAMETER :: PETSC_NORM_2 = NORM_2
  NormType, PARAMETER :: PETSC_NORM_INFINITY = NORM_INFINITY
  
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
#if ( PETSC_VERSION_GE(3,2,0) )
  PCType, PARAMETER ::  PETSC_PCGASM = PCGASM
  PCType, PARAMETER ::  PETSC_PCPARMS = PCPARMS
  PCType, PARAMETER ::  PETSC_PCTFS = PCTFS
#if ( PETSC_VERSION_LE(3,3,0) )
  PCType, PARAMETER ::  PETSC_PCPROMETHEUS = PCPROMETHEUS
#endif
  PCType, PARAMETER ::  PETSC_PCGALERKIN = PCGALERKIN
  PCType, PARAMETER ::  PETSC_PCEXOTIC = PCEXOTIC
#if ( PETSC_VERSION_LE(3,4,0) )
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
#if ( PETSC_VERSION_GE(3,2,0) )
#if ( PETSC_VERSION_LE(3,3,0) )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SPOOLES = MATSOLVERSPOOLES
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PLAPACK = MATSOLVERPLAPACK
#endif
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU = MATSOLVERSUPERLU
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU_DIST = MATSOLVERSUPERLU_DIST
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_UMFPACK = MATSOLVERUMFPACK
#if ( PETSC_VERSION_GE(3,5,0) )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_CHOLMOD = MATSOLVERCHOLMOD
#endif 
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_ESSL = MATSOLVERESSL
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_LUSOL = MATSOLVERLUSOL
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_MUMPS = MATSOLVERMUMPS
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PASTIX = MATSOLVERPASTIX
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_MATLAB = MATSOLVERMATLAB
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PETSC = MATSOLVERPETSC
#if ( PETSC_VERSION_GE(3,5,0) )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_BAS = MATSOLVERBAS
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_CUSPARSE = MATSOLVERCUSPARSE
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_BSTRM = MATSOLVERBSTRM
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_SBSTRM = MATSOLVERSBSTRM
#endif 
#else
#if ( PETSC_VERSION_GE(3,0,0) )
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
#if ( PETSC_VERSION_GE(3,1,0) )
  MatSolverPackage, PARAMETER :: PETSC_MAT_SOLVER_PASTIX = MAT_SOLVER_PASTIX
#endif
#endif
#endif
  
  !SNES converged types
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_FNORM_ABS = SNES_CONVERGED_FNORM_ABS
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_FNORM_RELATIVE = SNES_CONVERGED_FNORM_RELATIVE
#if ( PETSC_VERSION_LT(3,3,0) )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_PNORM_RELATIVE = SNES_CONVERGED_PNORM_RELATIVE
#endif
#if ( PETSC_VERSION_GE(3,5,0) )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_SNORM_RELATIVE = SNES_CONVERGED_SNORM_RELATIVE
#endif
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_ITS = SNES_CONVERGED_ITS
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_TR_DELTA = SNES_CONVERGED_TR_DELTA
#if ( PETSC_VERSION_GE(3,2,0) )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FUNCTION_DOMAIN = SNES_DIVERGED_FUNCTION_DOMAIN
#endif
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FUNCTION_COUNT = SNES_DIVERGED_FUNCTION_COUNT
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LINEAR_SOLVE = SNES_DIVERGED_LINEAR_SOLVE
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FNORM_NAN = SNES_DIVERGED_FNORM_NAN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_MAX_IT = SNES_DIVERGED_MAX_IT
#if ( PETSC_VERSION_GE(3,2,0) )
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LS_FAILURE = SNES_DIVERGED_LINE_SEARCH
#else
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LS_FAILURE = SNES_DIVERGED_LS_FAILURE
#endif
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LOCAL_MIN = SNES_DIVERGED_LOCAL_MIN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_ITERATING = SNES_CONVERGED_ITERATING
  
  !SNES types
#if ( PETSC_VERSION_GE(3,4,0) )
  SNESType, PARAMETER :: PETSC_SNESTEST = SNESTEST
  SNESType, PARAMETER :: PETSC_SNESLS = SNESNEWTONLS
  SNESType, PARAMETER :: PETSC_SNESTR = SNESNEWTONTR
#else
  SNESType, PARAMETER :: PETSC_SNESLS = SNESLS
  SNESType, PARAMETER :: PETSC_SNESTR = SNESTR
#endif
#if ( PETSC_VERSION_GE(3,2,0) )
  SNESType, PARAMETER :: PETSC_SNESPYTHON = SNESPYTHON
#endif
  SNESType, PARAMETER :: PETSC_SNESNRICHARDSON = SNESNRICHARDSON
  SNESType, PARAMETER :: PETSC_SNESKSPONLY = SNESKSPONLY
#if ( PETSC_VERSION_GE(3,3,0) )
  SNESType, PARAMETER :: PETSC_SNESVIRS = SNESVINEWTONRSLS
  SNESType, PARAMETER :: PETSC_SNESVISS = SNESVINEWTONSSLS
#endif
#if ( PETSC_VERSION_GE(3,3,0) )
  SNESType, PARAMETER :: PETSC_SNESNGMRES = SNESNGMRES
  SNESType, PARAMETER :: PETSC_SNESQN = SNESQN
  SNESType, PARAMETER :: PETSC_SNESSHELL = SNESSHELL
  SNESType, PARAMETER :: PETSC_SNESNCG = SNESNCG
  SNESType, PARAMETER :: PETSC_SNESFAS = SNESFAS
  SNESType, PARAMETER :: PETSC_SNESMS = SNESMS
#endif

#if ( PETSC_VERSION_GE(3,3,0) )
#if ( PETSC_VERSION_GE(3,5,0) )
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

#if ( PETSC_VERSION_GE(3,4,0) )
  !SNES QN types
  SNESQNType, PARAMETER :: PETSC_SNES_QN_LBFGS = SNES_QN_LBFGS 
  SNESQNType, PARAMETER :: PETSC_SNES_QN_BROYDEN = SNES_QN_BROYDEN
  SNESQNType, PARAMETER :: PETSC_SNES_QN_BADBROYDEN = SNES_QN_BADBROYDEN
  
  !SNES QN restart types
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_NONE = SNES_QN_RESTART_NONE 
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_POWELL = SNES_QN_RESTART_POWELL
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_PERIODIC = SNES_QN_RESTART_PERIODIC 

  !SNES QN scaling types
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_NONE = SNES_QN_SCALE_NONE 
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_SHANNO = SNES_QN_SCALE_SHANNO 
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_LINESEARCH = SNES_QN_SCALE_LINESEARCH  
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_JACOBIAN = SNES_QN_SCALE_JACOBIAN   
#endif

   !TS types
#if ( PETSC_VERSION_GE(3,2,0) )
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
#if ( PETSC_VERSION_LT(3,1,0) )
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
#if ( PETSC_VERSION_GE(3,2,0) )
  TSType, PARAMETER :: PETSC_TS_CRANK_NICHOLSON = TSCN
#else
  TSType, PARAMETER :: PETSC_TS_CRANK_NICHOLSON = TSCRANK_NICHOLSON
#endif  
  TSType, PARAMETER :: PETSC_TS_SUNDIALS = TSSUNDIALS
#if ( PETSC_VERSION_GE(3,2,0) )
  TSType, PARAMETER :: PETSC_TS_RUNGE_KUTTA = TSRK
#else
  TSType, PARAMETER :: PETSC_TS_RUNGE_KUTTA = TSRUNGE_KUTTA
#endif
#if ( PETSC_VERSION_GE(3,5,0) )
  TSType, PARAMETER :: PETSC_TS_PYTHON = TSPYTHON
#endif  
  TSType, PARAMETER :: PETSC_TS_THETA = TSTHETA
#if ( PETSC_VERSION_GE(3,5,0) )
  TSType, PARAMETER :: PETSC_TS_ALPHA = TSALPHA
#endif  
  TSType, PARAMETER :: PETSC_TS_GL = TSGL
#if ( PETSC_VERSION_GE(3,5,0) )
  TSType, PARAMETER :: PETSC_TS_SSP = TSSSP
  TSType, PARAMETER :: PETSC_TS_ARKIMEX = TSARKIMEX
  TSType, PARAMETER :: PETSC_TS_ROSW = TSROSW
  TSType, PARAMETER :: PETSC_TS_EIMEX = TSEIMEX
#endif  
#endif
#endif

  !TS convergence flags
#if ( PETSC_VERSION_GE(3,5,0) )
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
#if ( PETSC_VERSION_GE(3,0,0) )
  TSSundialsGramSchmidtType, PARAMETER :: PETSC_SUNDIALS_MODIFIED_GS = SUNDIALS_MODIFIED_GS
  TSSundialsGramSchmidtType, PARAMETER :: PETSC_SUNDIALS_CLASSICAL_GS = SUNDIALS_CLASSICAL_GS
#else
  TSSundialsGramSchmitdType, PARAMETER :: PETSC_SUNDIALS_MODIFIED_GS = SUNDIALS_MODIFIED_GS
  TSSundialsGramSchmitdType, PARAMETER :: PETSC_SUNDIALS_CLASSICAL_GS = SUNDIALS_CLASSICAL_GS
#endif
  
  !Module types

  !Module variables

  LOGICAL, SAVE :: petscHandleError

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
#if ( PETSC_VERSION_GE(3,2,0) )
      PetscBool flg
#else
      PetscTruth flg
#endif
      PetscInt ierr
    END SUBROUTINE KSPSetInitialGuessNonzero
    
#if ( PETSC_VERSION_GE(3,5,0) )
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
    
#if ( PETSC_VERSION_GE(3,5,0) )
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

    SUBROUTINE MatColoringApply(coloring,isColoring,ierr)
      MatColoring coloring
      ISColoring isColoring
      PetscInt ierr
    END SUBROUTINE MatColoringApply
    
    SUBROUTINE MatColoringCreate(A,coloring,ierr)
      Mat A
      MatColoring coloring
      PetscInt ierr
    END SUBROUTINE MatColoringCreate
    
    SUBROUTINE MatColoringDestroy(coloring,ierr)
      MatColoring coloring
      PetscInt ierr
    END SUBROUTINE MatColoringDestroy
    
    SUBROUTINE MatColoringSetFromOptions(coloring,ierr)
      MatColoring coloring
      PetscInt ierr
    END SUBROUTINE MatColoringSetFromOptions
    
    SUBROUTINE MatColoringSetType(coloring,coloringType,ierr)
      MatColoring coloring
      MatColoringType coloringType
      PetscInt ierr
    END SUBROUTINE MatColoringSetType
    
    SUBROUTINE MatCreate(comm,A,ierr)
      MPI_Comm comm
      Mat A
      PetscInt ierr
    END SUBROUTINE MatCreate

#if ( PETSC_VERSION_GE(3,3,0) )
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
        
#if ( PETSC_VERSION_GE(3,3,0) )
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

    SUBROUTINE MatDenseGetArrayF90(A,array,ierr)
      Mat A
      PetscScalar, POINTER :: array(:,:)
      PetscInt ierr
    END SUBROUTINE MatDenseGetArrayF90

    SUBROUTINE MatDenseRestoreArrayF90(A,array,ierr)
      Mat A
      PetscScalar, POINTER :: array(:,:)
      PetscInt ierr
    END SUBROUTINE MatDenseRestoreArrayF90

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

    SUBROUTINE MatFDColoringSetUp(A,iscoloring,fdcoloring,ierr)
      Mat A
      ISColoring iscoloring
      MatFDColoring fdcoloring
      PetscInt ierr
    END SUBROUTINE MatFDColoringSetUp

#if ( PETSC_VERSION_LT(3,6,0) )
    SUBROUTINE MatGetColoring(A,coloring_type,iscoloring,ierr)
      Mat A
      MatColoringType coloring_type
      ISColoring iscoloring 
      PetscInt ierr
    END SUBROUTINE MatGetColoring
#endif    
    
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
    
    SUBROUTINE MatRestoreRow(A,row,ncols,cols,values,ierr)
      Mat A
      PetscInt row
      PetscInt ncols
      PetscInt cols(*)
      PetscScalar values(*)
      PetscInt ierr
    END SUBROUTINE MatRestoreRow
    
    SUBROUTINE MatSeqAIJGetArrayF90(A,array,ierr)
      Mat A
      PetscScalar, POINTER :: array(:,:)
      PetscInt ierr
    END SUBROUTINE MatSeqAIJGetArrayF90

    SUBROUTINE MatSeqAIJGetMaxRowNonzeros(A,maxNumberNonZeros,ierr)
      Mat A
      PetscInt maxNumberNonZeros
      PetscInt ierr
    END SUBROUTINE MatSeqAIJGetMaxRowNonzeros

    SUBROUTINE MatSeqAIJRestoreArrayF90(A,array,ierr)
      Mat A
      PetscScalar, POINTER :: array(:,:)
      PetscInt ierr
    END SUBROUTINE MatSeqAIJRestoreArrayF90

    SUBROUTINE MatSetLocalToGlobalMapping(A,ctx,ierr)
      Mat A
      ISLocalToGlobalMapping ctx
      PetscInt ierr
    END SUBROUTINE MatSetLocalToGlobalMapping

    SUBROUTINE MatSetOption(A,option,flag,ierr)
      Mat A
      MatOption option
#if ( PETSC_VERSION_GE(3,2,0) )
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

#if ( PETSC_VERSION_GE(3,0,0) )
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

#if ( PETSC_VERSION_GE(3,4,0) )
    SUBROUTINE MatMumpsSetCntl(A,icntl,val,ierr)
      Mat A
      PetscInt icntl
      PetscReal val
      PetscInt ierr
    END SUBROUTINE MatMumpsSetCntl
#endif
#if ( PETSC_VERSION_MINOR >= 5 )
    SUBROUTINE PCSetReusePreconditioner(pc,flag,ierr)
      PC pc
      PetscBool flag
      PetscInt ierr
    END SUBROUTINE PCSetReusePreconditioner
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

#if ( PETSC_VERSION_GE(3,2,0) )
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

    SUBROUTINE SNESGetApplicationContext(snes,ctx,ierr)
      USE TYPES
      SNES snes
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE SNESGetApplicationContext

    SUBROUTINE SNESSetApplicationContext(snes,ctx,ierr)
      USE TYPES
      SNES snes
      TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE SNESSetApplicationContext

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
    
    SUBROUTINE SNESGetSolutionUpdate(snes,solutionUpdate,ierr)
      SNES snes
      Vec solutionUpdate
      PetscInt ierr
    END SUBROUTINE SNESGetSolutionUpdate

#if ( PETSC_VERSION_LT(3,6,0) )
    SUBROUTINE SNESSetFunctionNorm(snes,fnorm,ierr)
      SNES snes
      PetscReal fnorm
      PetscInt ierr
    END SUBROUTINE SNESSetFunctionNorm
#endif    

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
    
#if ( PETSC_VERSION_GE(3,2,0) )
    SUBROUTINE SnesLineSearchSetMonitor(linesearch,flag,ierr)
      SNESLineSearch linesearch
      PetscBool flag
      PetscInt ierr
    END SUBROUTINE SnesLineSearchSetMonitor
#endif

#if ( PETSC_VERSION_GE(3,3,0) )
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
    
#if ( PETSC_VERSION_LT(3,3,0) )
    SUBROUTINE SNESLineSearchSetParams(snes,alpha,maxstep,steptol,ierr)
      SNES snes
      PetscReal alpha
      PetscReal maxstep
      PetscReal steptol
      PetscInt ierr
    END SUBROUTINE SNESLineSearchSetParams
#endif    
    
#if ( PETSC_VERSION_GE(3,3,0) )
#if ( PETSC_VERSION_GE(3,4,0) )
    SUBROUTINE SNESGetLineSearch(snes,linesearch,ierr)
#else
    SUBROUTINE SNESGetSNESLineSearch(snes,linesearch,ierr)
#endif
      SNES snes
      SNESLineSearch linesearch
      PetscInt ierr
#if ( PETSC_VERSION_GE(3,4,0) )
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

    SUBROUTINE SNESGetFunction(snes,f,ffunction,ctx,ierr)
      USE TYPES
      SNES snes
      Vec f
      EXTERNAL ffunction
      PetscInt ctx
      !TYPE(SOLVER_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE SNESGetFunction

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

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
    SUBROUTINE SNESQNSetType(snes,qtype,ierr)
      SNES snes
      SNESQNType qtype
      PetscInt ierr
    END SUBROUTINE SNESQNSetType

    SUBROUTINE SNESQNSetRestartType(snes,rtype,ierr)
      SNES snes
      SNESQNRestartType rtype
      PetscInt ierr
    END SUBROUTINE SNESQNSetRestartType

    SUBROUTINE SNESQNSetScaleType(snes,stype,ierr)
      SNES snes
      SNESQNScaleType stype
      PetscInt ierr
    END SUBROUTINE SNESQNSetScaleType
#endif

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

#if ( PETSC_VERSION_LT(3,6,0) )
    SUBROUTINE SNESSetNormSchedule(snes,normschedule,ierr)
      SNES snes
      SNESNormSchedule normschedule
      PetscInt ierr
    END SUBROUTINE SNESSetNormSchedule
#else
    SUBROUTINE SNESSetNormType(snes,normtype,ierr)
      SNES snes
#if ( PETSC_VERSION_GE(3,5,0) )
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

    SUBROUTINE TSSetExactFinalTime(ts,eftopt,ierr)
      TS ts
      PetscBool eftopt
      PetscInt ierr
    END SUBROUTINE TSSetExactFinalTime

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
    
    SUBROUTINE TSSetRHSFunction(ts,r,rhsfunc,ctx,ierr)
      USE TYPES
      TS ts
      Vec r
      EXTERNAL rhsfunc
      TYPE(CELLML_PETSC_CONTEXT_TYPE), POINTER :: ctx
      PetscInt ierr
    END SUBROUTINE TSSetRHSFunction

    SUBROUTINE TSSetSolution(ts,initialsolution,ierr)
      TS ts
      Vec initialsolution
      PetscInt ierr
    END SUBROUTINE TSSetSolution

    SUBROUTINE TSGetSolution(ts,currentsolution,ierr)
      TS ts
      Vec currentsolution
      PetscInt ierr
    END SUBROUTINE TSGetSolution
    
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
    
    SUBROUTINE TSSolve(ts,x,ftime,ierr)
      TS ts
      Vec x
      PetscReal ftime
      PetscInt ierr
    END SUBROUTINE TSSolve

    SUBROUTINE TSStep(ts,steps,ptime,ierr)
      TS ts
      PetscInt steps
      PetscReal ptime
      PetscInt ierr
    END SUBROUTINE TSStep

    SUBROUTINE TSSundialsSetType(ts,sundialstype,ierr)
      TS ts
      TSSundialsType sundialstype
      PetscInt ierr
    END SUBROUTINE TSSundialsSetType

    SUBROUTINE TSSundialsSetTolerance(ts,abstol,reltol,ierr)
      TS ts
      PetscReal abstol
      PetscReal reltol
      PetscInt ierr
    END SUBROUTINE TSSundialsSetTolerance

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
      PetscInt ierr
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

    SUBROUTINE VecGetArrayReadF90(x,vec_data,ierr)
      Vec x
      PetscScalar, POINTER :: vec_data(:)
      PetscInt ierr
    END SUBROUTINE VecGetArrayReadF90

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

    SUBROUTINE VecNorm(x,ntype,val,ierr)
      Vec x
      NormType ntype
      PetscReal val
      PetscInt ierr
    END SUBROUTINE VecNorm

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

    SUBROUTINE VecRestoreArrayReadF90(x,vec_data,ierr)
      Vec x
      PetscScalar, POINTER :: vec_data(:)
      PetscInt ierr
    END SUBROUTINE VecRestoreArrayReadF90
    
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
#if ( PETSC_VERSION_LT(3,4,0) )
  PUBLIC PETSC_NULL
#endif
#if ( PETSC_VERSION_LT(3,3,0) )
  PUBLIC PETSC_NULL_TRUTH
#else
  PUBLIC PETSC_NULL_BOOL
#endif

  PUBLIC PETSC_ADD_VALUES,PETSC_INSERT_VALUES,PETSC_COMM_WORLD,PETSC_COMM_SELF,PETSC_DECIDE, &
    & PETSC_SCATTER_FORWARD,PETSC_SCATTER_REVERSE

#if ( PETSC_VERSION_GE(3,5,0) )
  PUBLIC PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_REAL
#else
  PUBLIC PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_DOUBLE_PRECISION
#endif  
 
  PUBLIC PETSC_KSPRICHARDSON,PETSC_KSPCG,PETSC_KSPCGNE,PETSC_KSPSTCG,PETSC_KSPGMRES,PETSC_KSPFGMRES, &
    & PETSC_KSPLGMRES,PETSC_KSPTCQMR,PETSC_KSPBCGS,PETSC_KSPBCGSL,PETSC_KSPCGS,PETSC_KSPTFQMR,PETSC_KSPCR,PETSC_KSPLSQR, &
    & PETSC_KSPPREONLY,PETSC_KSPQCG,PETSC_KSPBICG,PETSC_KSPMINRES,PETSC_KSPSYMMLQ,PETSC_KSPLCD

#if ( PETSC_VERSION_LT(3,3,0) )
  PUBLIC PETSC_KSPCHEBYCHEV
#else
  PUBLIC PETSC_KSPCHEBYSHEV
#endif

  PUBLIC PETSC_PCNONE,PETSC_PCJACOBI,PETSC_PCSOR,PETSC_PCLU,PETSC_PCSHELL,PETSC_PCBJACOBI,PETSC_PCMG,PETSC_PCEISENSTAT, &
    & PETSC_PCILU,PETSC_PCICC,PETSC_PCASM,PETSC_PCKSP,PETSC_PCCOMPOSITE,PETSC_PCREDUNDANT,PETSC_PCSPAI, &
    & PETSC_PCNN,PETSC_PCCHOLESKY,PETSC_PCPBJACOBI,PETSC_PCMAT,PETSC_PCHYPRE,PETSC_PCFIELDSPLIT,PETSC_PCML
#if ( PETSC_VERSION_GE(3,2,0) )
#if ( PETSC_VERSION_LT(3,4,0) )
  PUBLIC PETSC_PCPROMETHEUS
#endif
#if ( PETSC_VERSION_LE(3,4,0) )
  PUBLIC PETSC_PCHMPI,PETSC_PCASA
#endif  
  PUBLIC PETSC_PCGASM,PETSC_PCPARMS,PETSC_PCTFS,PETSC_PCGALERKIN,PETSC_PCEXOTIC, &
    & PETSC_PCSUPPORTGRAPH,PETSC_PCCP,PETSC_PCVFVT,PETSC_PCLSC,PETSC_PCPYTHON,PETSC_PCPFMG,PETSC_PCSYSPFMG, &
    & PETSC_PCREDISTRIBUTE,PETSC_PCSACUSP,PETSC_PCSACUSPPOLY,PETSC_PCBICGSTABCUSP,PETSC_PCSVD,PETSC_PCAINVCUSP
#else
  PUBLIC PETSC_PCMILU
#endif
  
#if ( PETSC_VERSION_LE(3,4,0) )
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
  
#if ( PETSC_VERSION_GE(3,5,0) )
  PUBLIC PETSc_KspSetReusePreconditioner
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
    & PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR
  
#if ( PETSC_VERSION_LT(3,4,0) )
  PUBLIC PETSC_NUM_MAT_OPTIONS
#endif
#if ( PETSC_VERSION_LT(3,5,0) )
  PUBLIC PETSC_MAT_CHECK_COMPRESSED_ROW
#endif
  PUBLIC PETSC_MAT_SPD,PETSC_MAT_NO_OFF_PROC_ENTRIES,PETSC_MAT_NO_OFF_PROC_ZERO_ROWS
  
  PUBLIC PETSC_MATCOLORING_NATURAL,PETSC_MATCOLORING_SL,PETSC_MATCOLORING_LF,PETSC_MATCOLORING_ID

  PUBLIC Petsc_MatColoringApply,Petsc_MatColoringCreate,Petsc_MatColoringDestroy,Petsc_MatColoringFinalise, &
    & Petsc_MatColoringInitialise,Petsc_MatColoringSetFromOptions,Petsc_MatColoringSetType
#if ( PETSC_VERSION_LT(3,6,0) )
  PUBLIC Petsc_MatGetColouring
#endif

  PUBLIC Petsc_MatFDColoringCreate,Petsc_MatFDColoringDestroy,Petsc_MatFDColoringFinalise,Petsc_MatFDColoringInitialise, &
    & Petsc_MatFDColoringSetFromOptions,Petsc_MatFDColoringSetFunction,Petsc_MatFDColoringSetParameters, &
    & Petsc_MatFDColoringSetup
#if ( PETSC_VERSION_LT(3,0,0) )
  PUBLIC Petsc_MatFDColoringSetFunctionSnes
#endif

  PUBLIC PETSC_MATINITIALISE,PETSC_MATFINALISE,PETSC_MATASSEMBLYBEGIN,PETSC_MATASSEMBLYEND,PETSC_MATCREATE, &
    & PETSC_MATCREATESEQAIJ,PETSC_MATCREATESEQDENSE,PETSC_MATDESTROY, &
    & PETSC_MATGETINFO,PETSC_MATGETOWNERSHIPRANGE, &
    & PETSC_MATGETROW,PETSC_MATGETVALUES,PETSC_MATRESTOREROW, &
    & PETSC_MATSETLOCALTOGLOBALMAPPING,PETSC_MATSETOPTION,PETSC_MATSETSIZES,PETSC_MATSETVALUE,PETSC_MATSETVALUES, &
    & PETSC_MATSETVALUELOCAL,PETSC_MATSETVALUESLOCAL,PETSC_MATVIEW,PETSC_MATZEROENTRIES,PETSC_MATSETTYPE
#if ( PETSC_VERSION_GE(3,3,0) )
  PUBLIC Petsc_MatDenseGetArrayF90,Petsc_MatDenseRestoreArrayF90
  
  PUBLIC Petsc_MatSeqAIJGetArrayF90,Petsc_MatSeqAIJGetMaxRowNonzeros,Petsc_MatSeqAIJRestoreArrayF90
#endif
#if ( PETSC_VERSION_GE(3,3,0) )
  PUBLIC PETSC_MATCREATEAIJ,PETSC_MATCREATEDENSE
#else
  PUBLIC PETSC_MATCREATEMPIAIJ,PETSC_MATCREATEMPIDENSE
#endif
  
#if ( PETSC_VERSION_GE(3,0,0) )
#if ( PETSC_VERSION_LT(3,4,0) )
  PUBLIC PETSC_MAT_SOLVER_SPOOLES,PETSC_MAT_SOLVER_PLAPACK
#endif
  PUBLIC PETSC_MAT_SOLVER_SUPERLU,PETSC_MAT_SOLVER_SUPERLU_DIST,PETSC_MAT_SOLVER_UMFPACK, &
    & PETSC_MAT_SOLVER_ESSL,PETSC_MAT_SOLVER_LUSOL,PETSC_MAT_SOLVER_MUMPS,PETSC_MAT_SOLVER_MATLAB, &
    & PETSC_MAT_SOLVER_PETSC
#if ( PETSC_VERSION_GE(3,1,0) )
  PUBLIC PETSC_MAT_SOLVER_PASTIX
#endif
#endif
    
  PUBLIC PETSC_NORM_1,PETSC_NORM_2,PETSC_NORM_INFINITY

  PUBLIC PETSC_PCINITIALISE,PETSC_PCFINALISE,PETSC_PCSETTYPE

#if ( PETSC_VERSION_GE(3,0,0) )
  PUBLIC PETSC_PCFACTORSETMATSOLVERPACKAGE
  PUBLIC Petsc_PCFactorSetUpMatSolverPackage
  PUBLIC Petsc_PCFactorGetMatrix
  PUBLIC Petsc_MatMumpsSetIcntl
#if ( PETSC_VERSION_GE(3,4,0) )
  PUBLIC Petsc_MatMumpsSetCntl
#endif
#if ( PETSC_VERSION_MINOR >= 5 )
  PUBLIC PETSC_PCSETREUSEPRECONDITIONER
#endif
#endif
  
  PUBLIC PETSC_TS_EULER,PETSC_TS_BEULER,PETSC_TS_PSEUDO,PETSC_TS_SUNDIALS,PETSC_TS_CRANK_NICHOLSON,PETSC_TS_RUNGE_KUTTA
#if ( PETSC_VERSION_GE(3,2,0) )
  PUBLIC PETSC_TS_PYTHON,PETSC_TS_THETA,PETSC_TS_ALPHA,PETSC_TS_SSP,PETSC_TS_ARKIMEX
#else
#if ( PETSC_VERSION_GE(3,1,0) )
  PUBLIC PETSC_TS_THETA,PETSC_TS_GL
#endif
#endif

  PUBLIC PETSC_TS_LINEAR,PETSC_TS_NONLINEAR

  PUBLIC PETSC_SUNDIALS_ADAMS,PETSC_SUNDIALS_BDF,PETSC_SUNDIALS_MODIFIED_GS,PETSC_SUNDIALS_CLASSICAL_GS

  PUBLIC PETSC_TSCREATE,PETSC_TSDESTROY,PETSC_TSFINALISE,PETSC_TSINITIALISE,PETSC_TSMONITORSET, &
    & PETSC_TSSETDURATION,PETSC_TSSETFROMOPTIONS, PETSC_TSSETEXACTFINALTIME, &
    & PETSC_TSSETINITIALTIMESTEP,PETSC_TSSETSOLUTION, PETSC_TSGETSOLUTION, &
    & PETSC_TSSETPROBLEMTYPE,PETSC_TSSETRHSFUNCTION,PETSC_TSSETTIMESTEP,PETSC_TSSETTYPE, &
    & PETSC_TSSUNDIALSSETTYPE, PETSC_TSSUNDIALSSETTOLERANCE,PETSC_TSSOLVE,PETSC_TSSTEP
#if ( PETSC_VERSION_LT(3,2,0) )
  PUBLIC PETSC_TSSETMATRICES
#endif
  
  PUBLIC PETSC_ERRORHANDLING_SET_OFF,PETSC_ERRORHANDLING_SET_ON
  
  PUBLIC PETSC_FINALIZE,PETSC_INITIALIZE
#if ( PETSC_VERSION_GE(3,2,0) )
  PUBLIC PETSC_LOGVIEW
#else
  PUBLIC PETSC_LOGPRINTSUMMARY
#endif
  
  PUBLIC PETSC_SNESLS,PETSC_SNESTR,PETSC_SNESTEST
#if ( PETSC_VERSION_GE(3,2,0) )
  PUBLIC PETSC_SNESPYTHON
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 4 )
  PUBLIC PETSC_SNESQN
#endif
  
#if ( PETSC_VERSION_GE(3,3,0) )
  PUBLIC PETSC_SNES_LINESEARCH_BASIC,PETSC_SNES_LINESEARCH_BT,PETSC_SNES_LINESEARCH_L2,PETSC_SNES_LINESEARCH_CP, &
    & PETSC_SNES_LINESEARCH_SHELL
  
  PUBLIC PETSC_SNES_LINESEARCH_LINEAR,PETSC_SNES_LINESEARCH_QUADRATIC,PETSC_SNES_LINESEARCH_CUBIC
#else
  PUBLIC PETSC_SNES_LINESEARCH_NONORMS,PETSC_SNES_LINESEARCH_NO,PETSC_SNES_LINESEARCH_QUADRATIC,PETSC_SNES_LINESEARCH_CUBIC
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 3 )

  PUBLIC PETSC_SNES_NORM_DEFAULT,PETSC_SNES_NORM_NONE,PETSC_SNES_NORM_INITIAL_ONLY, &
    & PETSC_SNES_NORM_FINAL_ONLY,PETSC_SNES_NORM_INITIAL_FINAL_ONLY
#if ( PETSC_VERSION_GE(3,5,0) )
  PUBLIC PETSC_SNES_NORM_ALWAYS
#else
  PUBLIC PETSC_SNES_NORM_FUNCTION
#endif
#endif

  PUBLIC PETSC_SNES_CONVERGED_FNORM_ABS,PETSC_SNES_CONVERGED_FNORM_RELATIVE, &
    & PETSC_SNES_CONVERGED_ITS,PETSC_SNES_CONVERGED_TR_DELTA,PETSC_SNES_DIVERGED_FUNCTION_COUNT,PETSC_SNES_DIVERGED_LINEAR_SOLVE, &
    & PETSC_SNES_DIVERGED_FNORM_NAN,PETSC_SNES_DIVERGED_MAX_IT,PETSC_SNES_DIVERGED_LOCAL_MIN,PETSC_SNES_CONVERGED_ITERATING
#if ( PETSC_VERSION_GE(3,2,0) )
  PUBLIC PETSC_SNES_DIVERGED_FUNCTION_DOMAIN,PETSC_SNES_DIVERGED_LS_FAILURE
#else
  PUBLIC PETSC_SNES_DIVERGED_LS_FAILURE
#endif
#if ( PETSC_VERSION_LE(3,2,0) )
  PUBLIC PETSC_SNES_CONVERGED_PNORM_RELATIVE
#endif

#if ( PETSC_VERSION_GT(3,6,0) )
  PUBLIC Petsc_SnesGetApplicationContext,Petsc_SnesSetApplicationContext
#endif
  
  PUBLIC PETSC_SNESFINALISE,PETSC_SNESINITIALISE,PETSC_SNESCREATE,PETSC_SNESDESTROY,PETSC_SNESGETCONVERGEDREASON, &
    & petsc_sneslinesearchsetnorms,petsc_sneslinesearchgetnorms, &
    & PETSC_SNESGETITERATIONNUMBER,PETSC_SNESGETKSP,PETSC_SNESMONITORSET,PETSC_SNESSETFROMOPTIONS,PETSC_SNESSETFUNCTION, &
    & PETSC_SNESSETJACOBIAN,PETSC_SNESSETTOLERANCES,PETSC_SNESSETTRUSTREGIONTOLERANCE,PETSC_SNESSETTYPE,PETSC_SNESSOLVE, &
    & PETSC_SNESSETKSP,PETSC_SNESGETJACOBIAN,Petsc_SnesComputeJacobianDefaultColor,Petsc_SnesComputeJacobianDefault, &
    & PETSC_SNESSETCONVERGENCETEST,Petsc_SnesLineSearchGetVecs,Petsc_SnesGetFunction,Petsc_SnesGetSolutionUpdate
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  PUBLIC PETSC_SNESSETNORMSCHEDULE
#else
  PUBLIC PETSC_SNESSETNORMTYPE
#endif
#if ( PETSC_VERSION_GE(3,2,0) )
  PUBLIC Petsc_SnesLineSearchSetMonitor
#endif
#if ( PETSC_VERSION_LT(3,6,0) )
  PUBLIC Petsc_SnesSetFunctionNorm,Petsc_SnesSetNormType
#endif  
#if ( PETSC_VERSION_GE(3,3,0) )
  PUBLIC Petsc_SnesLineSearchFinalise,Petsc_SnesLineSearchInitialise
  PUBLIC Petsc_SnesLineSearchSetComputeNorms,Petsc_SnesLineSearchComputeNorms, &
    & Petsc_SnesLineSearchSetOrder,Petsc_SnesLineSearchSetType, &
    & Petsc_SnesLineSearchBTSetAlpha,Petsc_SnesLineSearchSetTolerances
#else
  PUBLIC PETSC_SNESLINESEARCHSET,PETSC_SNESLINESEARCHSETPARAMS
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR == 3 )
  PUBLIC Petsc_SnesGetSnesLineSearch
#elif ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 4 )
  PUBLIC Petsc_SnesGetLineSearch
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  PUBLIC PETSC_SNESQNSETTYPE,PETSC_SNESQNSETRESTARTTYPE,PETSC_SNESQNSETSCALETYPE
#endif
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 4 )
  PUBLIC PETSC_SNES_QN_LBFGS,PETSC_SNES_QN_BROYDEN,PETSC_SNES_QN_BADBROYDEN
  PUBLIC PETSC_SNES_QN_RESTART_NONE,PETSC_SNES_QN_RESTART_POWELL,PETSC_SNES_QN_RESTART_PERIODIC
  PUBLIC PETSC_SNES_QN_SCALE_NONE,PETSC_SNES_QN_SCALE_SHANNO,PETSC_SNES_QN_SCALE_LINESEARCH,PETSC_SNES_QN_SCALE_JACOBIAN
#endif
  
  PUBLIC PETSC_VECINITIALISE,PETSC_VECFINALISE,PETSC_VECASSEMBLYBEGIN,PETSC_VECASSEMBLYEND,PETSC_VECCOPY,PETSC_VECCREATE, &
    & PETSC_VECCREATEGHOST,PETSC_VECCREATEGHOSTWITHARRAY,PETSC_VECCREATEMPI,PETSC_VECCREATEMPIWITHARRAY,PETSC_VECCREATESEQ, &
    & PETSC_VECCREATESEQWITHARRAY,PETSC_VECDESTROY,PETSC_VECDUPLICATE,PETSC_VECGETARRAYF90, &
    & PETSC_VECGETLOCALSIZE,PETSC_VECGETOWNERSHIPRANGE,PETSC_VECGETSIZE,PETSC_VECGETVALUES,PETSC_VECGHOSTGETLOCALFORM, &
    & PETSC_VECGHOSTRESTORELOCALFORM,PETSC_VECGHOSTUPDATEBEGIN,PETSC_VECGHOSTUPDATEEND, &
    & PETSC_VECRESTOREARRAYF90,PETSC_VECSCALE,PETSC_VECSET,PETSC_VECSETFROMOPTIONS,PETSC_VECSETLOCALTOGLOBALMAPPING, &
    & PETSC_VECSETSIZES,PETSC_VECSETVALUES,PETSC_VECSETVALUESLOCAL,PETSC_VECVIEW,Petsc_VecDot,Petsc_VecNorm
#if ( PETSC_VERSION_GE(3,6,0) )
  PUBLIC Petsc_VecGetArrayReadF90,Petsc_VecRestoreArrayReadF90
#endif

  PUBLIC PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_DRAW_WORLD,PETSC_VIEWER_DRAW_SELF

CONTAINS

  !
  !================================================================================================================================
  !

  !>Set PETSc error handling on
  SUBROUTINE PETSC_ERRORHANDLING_SET_OFF(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ERRORHANDLING_SET_OFF",err,error,*999)

    petscHandleError=.FALSE.
    
    EXITS("PETSC_ERRORHANDLING_SET_OFF")
    RETURN
999 ERRORSEXITS("PETSC_ERRORHANDLING_SET_OFF",err,error)
    RETURN 1
  END SUBROUTINE PETSC_ERRORHANDLING_SET_OFF
    
  !
  !================================================================================================================================
  !

  !>Set PETSc error handling on
  SUBROUTINE PETSC_ERRORHANDLING_SET_ON(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ERRORHANDLING_SET_ON",err,error,*999)

    petscHandleError=.TRUE.
    
    EXITS("PETSC_ERRORHANDLING_SET_ON")
    RETURN
999 ERRORSEXITS("PETSC_ERRORHANDLING_SET_ON",err,error)
    RETURN 1
  END SUBROUTINE PETSC_ERRORHANDLING_SET_ON
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscFinalize routine
  SUBROUTINE PETSC_FINALIZE(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_FINALIZE",err,error,*999)

    CALL PetscFinalize(ERR)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscFinalize",err,error,*999)
    ENDIF
    
    EXITS("PETSC_FINALIZE")
    RETURN
999 ERRORSEXITS("PETSC_FINALIZE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_FINALIZE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscInitialize routine.
  SUBROUTINE PETSC_INITIALIZE(FILE,err,error,*)

    !Argument Variables
    CHARACTER(LEN=*), INTENT(IN) :: FILE !<Filename for PETSc options file
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_INITIALIZE",err,error,*999)

    petscHandleError=.TRUE.
    CALL PetscInitialize(FILE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscInitialize",err,error,*999)
    ENDIF
    ! Disable PETSc's signal handler as we have our own OpenCMISS signal handlers in cmiss_c.c
    CALL PetscPopSignalHandler(ERR)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscPopSignalHandler",err,error,*999)
    ENDIF
    
    EXITS("PETSC_INITIALIZE")
    RETURN
999 ERRORSEXITS("PETSC_INITIALIZE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_INITIALIZE
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc IS structure and destroy the IS
  SUBROUTINE PETSC_ISFINALISE(is,err,error,*)

    !Argument Variables
    TYPE(PetscISType), INTENT(INOUT) :: is !<The IS to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISFINALISE",err,error,*999)

    IF(is%is/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_ISDESTROY(is,err,error,*999)
    ENDIF
    
    EXITS("PETSC_ISFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_ISFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc IS structure
  SUBROUTINE PETSC_ISINITIALISE(is,err,error,*)

    !Argument Variables
    TYPE(PetscISType), INTENT(INOUT) :: is !<The IS to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISINITIALISE",err,error,*999)

    is%is=PETSC_NULL_OBJECT
    
    EXITS("PETSC_ISINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_ISINITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISDestroy routine.
  SUBROUTINE PETSC_ISDESTROY(is,err,error,*)

    !Argument Variables
    TYPE(PetscISType), INTENT(INOUT) :: is !<The index set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISDESTROY",err,error,*999)

    CALL ISDestroy(is%is,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISDestroy",err,error,*999)
    ENDIF
    is%is=PETSC_NULL_OBJECT
        
    EXITS("PETSC_ISDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_ISDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISDESTROY
    
  !
  !
  !================================================================================================================================
  !

  !Finalise the PETSc ISColoring structure and destroy the ISColoring
  SUBROUTINE PETSC_ISCOLORINGFINALISE(iscoloring,err,error,*)

    !Argument Variables
    TYPE(PetscISColoringType), INTENT(INOUT) :: iscoloring !<The ISColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISCOLORINGFINALISE",err,error,*999)

    IF(iscoloring%iscoloring/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_ISCOLORINGDESTROY(iscoloring,err,error,*999)
    ENDIF
    
    EXITS("PETSC_ISCOLORINGFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_ISCOLORINGFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISCOLORINGFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc ISColoring structure
  SUBROUTINE Petsc_ISColoringInitialise(iscoloring,err,error,*)

    !Argument Variables
    TYPE(PetscISColoringType), INTENT(INOUT) :: iscoloring !<The ISColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISColoringInitialise",err,error,*999)

    iscoloring%iscoloring=PETSC_NULL_OBJECT
    
    EXITS("Petsc_ISColoringInitialise")
    RETURN
999 ERRORSEXITS("Petsc_ISColoringInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISColoringInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISColoringDestroy routine.
  SUBROUTINE PETSC_ISCOLORINGDESTROY(iscoloring,err,error,*)

    !Argument Variables
    TYPE(PetscISColoringType), INTENT(INOUT) :: iscoloring !<The index set coloring
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISCOLORINGDESTROY",err,error,*999)

    CALL ISColoringDestroy(iscoloring%iscoloring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISColoringDestroy",err,error,*999)
    ENDIF
    iscoloring%iscoloring=PETSC_NULL_OBJECT
    
    EXITS("PETSC_ISCOLORINGDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_ISCOLORINGDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISCOLORINGDESTROY
  
  !  
  !================================================================================================================================
  !

  !Finalise the PETSc ISLocalToGlobalMapping structure and destroy the ISLocalToGlobalMapping
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGFINALISE(ISLOCALTOGLOBALMAPPING,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: ISLOCALTOGLOBALMAPPING !<The ISLocalToGlobalMapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE",err,error,*999)

    IF(ISLOCALTOGLOBALMAPPING%ISLOCALTOGLOBALMAPPING/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_ISLOCALTOGLOBALMAPPINGDESTROY(ISLOCALTOGLOBALMAPPING,err,error,*999)
    ENDIF
    
    EXITS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_ISLOCALTOGLOBALMAPPINGFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGFINALISE
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc ISLocalToGlobalMapping structure
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE(ISLOCALTOGLOBALMAPPING,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: ISLOCALTOGLOBALMAPPING !<The ISLocalToGlobalMapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE",err,error,*999)

    ISLOCALTOGLOBALMAPPING%ISLOCALTOGLOBALMAPPING=PETSC_NULL_OBJECT
    
    EXITS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingApply routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLY(CTX,TYPE,NIN,IDXIN,NOUT,IDXOUT,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(IN) :: TYPE !<The type of local to global mapping
    INTEGER(INTG), INTENT(IN) :: NIN !<The number of local indicies
    INTEGER(INTG), INTENT(IN) :: IDXIN(*)
    INTEGER(INTG), INTENT(OUT) :: NOUT
    INTEGER(INTG), INTENT(OUT) :: IDXOUT(*)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY",err,error,*999)

    CALL ISLocalToGlobalMappingApply(CTX%ISLOCALTOGLOBALMAPPING,TYPE,NIN,IDXIN,NOUT,IDXOUT,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingApply",err,error,*999)
    ENDIF
    CTX%ISLOCALTOGLOBALMAPPING=PETSC_NULL_OBJECT
    
    EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY")
    RETURN
999 ERRORSEXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingApplyIS routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS(CTX,ISIN,ISOUT,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: CTX !<The local to global mapping context
    TYPE(PetscISType), INTENT(IN) :: ISIN
    TYPE(PetscISType), INTENT(OUT) :: ISOUT
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS",err,error,*999)

    CALL ISLocalToGlobalMappingApplyIS(CTX%ISLOCALTOGLOBALMAPPING,ISIN%is,ISOUT%is,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingApplyIS",err,error,*999)
    ENDIF
    
    EXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS")
    RETURN
999 ERRORSEXITS("PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS",err,error)
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGAPPLYIS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingCreate routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGCREATE(COMMUNICATOR,N,GLOBALNUM,CTX,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: N !<The number of local indices
    INTEGER(INTG), INTENT(IN) :: GLOBALNUM(*) !<The global number for each local index
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE",err,error,*999)

    CALL ISLocalToGlobalMappingCreate(COMMUNICATOR,N,GLOBALNUM,CTX%ISLOCALTOGLOBALMAPPING,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingCreate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE")
    RETURN
999 ERRORSEXITS("PETSC_ISLOCALTOGLOBALMAPPINGCREATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingDestroy routine.
  SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGDESTROY(CTX,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY",err,error,*999)

    CALL ISLocalToGlobalMappingDestroy(CTX%ISLOCALTOGLOBALMAPPING,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingDestroy",err,error,*999)
    ENDIF
    CTX%ISLOCALTOGLOBALMAPPING=PETSC_NULL_OBJECT
    
    EXITS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_ISLOCALTOGLOBALMAPPINGDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_ISLOCALTOGLOBALMAPPINGDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPCreate routine
  SUBROUTINE PETSC_KSPCREATE(COMMUNICATOR,ksp,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator for the KSP creation
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<On exit, the Ksp information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPCREATE",err,error,*999)

    CALL KSPCreate(COMMUNICATOR,ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPCreate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPCREATE")
    RETURN
999 ERRORSEXITS("PETSC_KSPCREATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPDestroy routine
  SUBROUTINE PETSC_KSPDESTROY(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPDESTROY",err,error,*999)

    CALL KSPDestroy(ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPDestroy",err,error,*999)
    ENDIF
    ksp%ksp=PETSC_NULL_OBJECT
    
    EXITS("PETSC_KSPDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_KSPDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_KSPDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetConvergedReason routine
  SUBROUTINE PETSC_KSPGETCONVERGEDREASON(ksp,REASON,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP information
    INTEGER(INTG), INTENT(OUT) :: REASON !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPGETCONVERGEDREASON",err,error,*999)

    CALL KSPGetConvergedReason(ksp%ksp,REASON,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetConvergedReason",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPGETCONVERGEDREASON")
    RETURN
999 ERRORSEXITS("PETSC_KSPGETCONVERGEDREASON",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPGETCONVERGEDREASON
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetIterationNumber routine
  SUBROUTINE PETSC_KSPGETITERATIONNUMBER(ksp,ITERATION_NUMBER,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP information
    INTEGER(INTG), INTENT(OUT) :: ITERATION_NUMBER !<On exit, the number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPGETITERATIONNUMBER",err,error,*999)

    CALL KSPGetIterationNumber(ksp%ksp,ITERATION_NUMBER,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetIterationNumber",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPGETITERATIONNUMBER")
    RETURN
999 ERRORSEXITS("PETSC_KSPGETITERATIONNUMBER",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPGETITERATIONNUMBER
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetPC routine
  SUBROUTINE PETSC_KSPGETPC(ksp,pc,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to get the PC for
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<On exit, the PC associated with the Ksp
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPGETPC",err,error,*999)

    CALL KSPGetPC(ksp%ksp,pc%pc,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetPC",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPGETPC")
    RETURN
999 ERRORSEXITS("PETSC_KSPGETPC",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPGETPC
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetResidualNorm routine
  SUBROUTINE PETSC_KSPGETRESIDUALNORM(ksp,RESIDUAL_NORM,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to get the PC for
    REAL(DP), INTENT(OUT) :: RESIDUAL_NORM !<On exit, the residual norm for the KSP
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPGETRESIDUALNORM",err,error,*999)

    CALL KSPGetResidualNorm(ksp%ksp,RESIDUAL_NORM,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetResidualNorm.",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPGETRESIDUALNORM")
    RETURN
999 ERRORSEXITS("PETSC_KSPGETRESIDUALNORM",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPGETRESIDUALNORM
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc KSP structure and destroy the KSP
  SUBROUTINE PETSC_KSPFINALISE(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPFINALISE",err,error,*999)

    IF(ksp%ksp/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_KSPDESTROY(ksp,err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_KSPFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_KSPFINALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGMRESSetRestart routine
  SUBROUTINE PETSC_KSPGMRESSETRESTART(ksp,RESTART,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the GMRES restart for
    INTEGER(INTG), INTENT(OUT) :: RESTART !<The restart value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPGMRESSETRESTART",err,error,*999)

    CALL KSPGMRESSetRestart(ksp%ksp,RESTART,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGMRESSetRestart.",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPGMRESSETRESTART")
    RETURN
999 ERRORSEXITS("PETSC_KSPGMRESSETRESTART",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPGMRESSETRESTART
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc KSP structure
  SUBROUTINE PETSC_KSPINITIALISE(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPINITIALISE",err,error,*999)

    ksp%ksp=PETSC_NULL_OBJECT
    
    EXITS("PETSC_KSPINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_KSPINITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_KSPINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetFromOptions routine
  SUBROUTINE PETSC_KSPSETFROMOPTIONS(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSETFROMOPTIONS",err,error,*999)

    CALL KSPSetFromOptions(ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetFromOptions",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSETFROMOPTIONS")
    RETURN
999 ERRORSEXITS("PETSC_KSPSETFROMOPTIONS",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetInitialGuessNonzero routine
  SUBROUTINE PETSC_KSPSETINITIALGUESSNONZERO(ksp,FLAG,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the initial guess non zero for
    LOGICAL, INTENT(IN) :: FLAG
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSETINITIALGUESSNONZERO",err,error,*999)

    IF(FLAG) THEN
      CALL KSPSetInitialGuessNonzero(ksp%ksp,PETSC_TRUE,err)
    ELSE
      CALL KSPSetInitialGuessNonzero(ksp%ksp,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetInitialGuessNonzero",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSETINITIALGUESSNONZERO")
    RETURN
999 ERRORSEXITS("PETSC_KSPSETINITIALGUESSNONZERO",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPSETINITIALGUESSNONZERO
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,5,0) )
  !>Buffer routine to the PETSc KSPSetOperators routine
  SUBROUTINE PETSc_KSPSetOperators(ksp,amat,pmat,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the operators for
    TYPE(PetscMatType), INTENT(INOUT) :: amat !<The matrix associated with the linear system
    TYPE(PetscMatType), INTENT(INOUT) :: pmat !<The matrix to be used in constructing the preconditioner
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSetOperators",err,error,*999)

    CALL KSPSetOperators(ksp%ksp,amat%mat,pmat%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetFromOperators.",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSetOperators")
    RETURN
999 ERRORSEXITS("PETSC_KSPSetOperators",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_KSPSetOperators
  
#else
  
  !>Buffer routine to the PETSc KSPSetOperators routine
  SUBROUTINE PETSC_KSPSETOPERATORS(ksp,AMAT,PMAT,FLAG,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the operators for
    TYPE(PetscMatType), INTENT(INOUT) :: AMAT !<The matrix associated with the linear system
    TYPE(PetscMatType), INTENT(INOUT) :: PMAT !<The matrix to be used in constructing the preconditioner
    MatStructure, INTENT(IN) :: FLAG !<Preconditioner matrix structure flag
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSETOPERATORS",err,error,*999)

    CALL KSPSetOperators(ksp%ksp,AMAT%mat,PMAT%mat,FLAG,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetFromOperators",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSETOPERATORS")
    RETURN
999 ERRORSEXITS("PETSC_KSPSETOPERATORS",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_KSPSETOPERATORS
#endif
  
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,5,0) )
  !>Buffer routine to the PETSc KSPSetReusePreconditioner routine
  SUBROUTINE PETSc_KSPSetReusePreconditioner(ksp,flag,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the options for
    LOGICAL, INTENT(IN) :: flag
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSc_KSPSetReusePreconditioner",err,error,*999)

    IF(flag) THEN
      CALL KSPSetReusePreconditioner(ksp%ksp,PETSC_TRUE,err)
    ELSE
      CALL KSPSetReusePreconditioner(ksp%ksp,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetReusePreconditioner.",err,error,*999)
    ENDIF
    
    EXITS("PETSc_KSPSetReusePreconditioner")
    RETURN
999 ERRORSEXITS("PETSc_KSPSetReusePreconditioner",err,error)
    RETURN 1
    
  END SUBROUTINE PETSc_KSPSetReusePreconditioner
#endif
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetTolerances routine
  SUBROUTINE PETSC_KSPSETTOLERANCES(ksp,RTOL,ATOL,DTOL,MAXITS,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the tolerances for
    REAL(DP), INTENT(IN) :: RTOL !<The relative tolerance to set
    REAL(DP), INTENT(IN) :: ATOL !<The absolution tolerance to set
    REAL(DP), INTENT(IN) :: DTOL !<The divergence tolerance to set
    INTEGER(INTG), INTENT(IN) :: MAXITS !<The maximum number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSETTOLERANCES",err,error,*999)

    CALL KSPSetTolerances(ksp%ksp,RTOL,ATOL,DTOL,MAXITS,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetTolerances",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSETTOLERANCES")
    RETURN
999 ERRORSEXITS("PETSC_KSPSETTOLERANCES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPSETTOLERANCES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetType routine
  SUBROUTINE PETSC_KSPSETTYPE(ksp,METHOD,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the type for
    KSPType, INTENT(IN) :: METHOD !<The Ksp method
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSETTYPE",err,error,*999)

    CALL KSPSetType(ksp%ksp,METHOD,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_KSPSETTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPSETTYPE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetUp routine
  SUBROUTINE PETSC_KSPSETUP(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set up
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSETUP",err,error,*999)

    CALL KSPSetUp(ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetUp",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSETUP")
    RETURN
999 ERRORSEXITS("PETSC_KSPSETUP",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPSETUP
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSolve routine
  SUBROUTINE PETSC_KSPSOLVE(ksp,B,X,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set up
    TYPE(PetscVecType), INTENT(INOUT) :: B !<The RHS vector
    TYPE(PetscVecType), INTENT(INOUT)  :: X !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_KSPSOLVE",err,error,*999)

    CALL KSPSolve(ksp%ksp,B%vec,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSolve",err,error,*999)
    ENDIF
    
    EXITS("PETSC_KSPSOLVE")
    RETURN
999 ERRORSEXITS("PETSC_KSPSOLVE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_KSPSOLVE
    
#if ( PETSC_VERSION_GE(3,2,0) )
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscLogView routine.
  SUBROUTINE PETSC_LOGVIEW(VIEWER,err,error,*)

    !Argument Variables
    PetscViewer, INTENT(IN) :: VIEWER !<The viewer to print the log to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_LOGVIEW",err,error,*999)

    CALL PetscLogView(VIEWER,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscLogView.",err,error,*999)
    ENDIF
    
    EXITS("PETSC_LOGVIEW")
    RETURN
999 ERRORSEXITS("PETSC_LOGVIEW",err,error)
    RETURN 1
  END SUBROUTINE PETSC_LOGVIEW
  
#else
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscLogPrintSummary routine.
  SUBROUTINE PETSC_LOGPRINTSUMMARY(COMMUNICATOR,FILE,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    CHARACTER(LEN=*), INTENT(IN) :: FILE !<Filename for the log summary
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_LOGPRINTSUMMARY",err,error,*999)

    CALL PetscLogPrintSummary(COMMUNICATOR,FILE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscLogPrintSummary",err,error,*999)
    ENDIF
    
    EXITS("PETSC_LOGPRINTSUMMARY")
    RETURN
999 ERRORSEXITS("PETSC_LOGPRINTSUMMARY",err,error)
    RETURN 1
  END SUBROUTINE PETSC_LOGPRINTSUMMARY
    
  !
  !================================================================================================================================
  !
  
#endif
  
  !>Buffer routine to the PETSc MatAssemblyBegin routine.
  SUBROUTINE PETSC_MATASSEMBLYBEGIN(A,ASSEMBLY_TYPE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !The matrix to assemble
    MatAssemblyType, INTENT(IN) :: ASSEMBLY_TYPE !<The assembly type 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATASSEMBLYBEGIN",err,error,*999)

    CALL MatAssemblyBegin(A%mat,ASSEMBLY_TYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatAssemblyBegin",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATASSEMBLYBEGIN")
    RETURN
999 ERRORSEXITS("PETSC_MATASSEMBLYBEGIN",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATASSEMBLYBEGIN
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatAssemblyEnd routine.
  SUBROUTINE PETSC_MATASSEMBLYEND(A,ASSEMBLY_TYPE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to assemble
    MatAssemblyType, INTENT(IN) :: ASSEMBLY_TYPE !<The assembly type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATASSEMBLYEND",err,error,*999)

    CALL MatAssemblyEnd(A%mat,ASSEMBLY_TYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatAssemblyEnd",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATASSEMBLYEND")
    RETURN
999 ERRORSEXITS("PETSC_MATASSEMBLYEND",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATASSEMBLYEND
       
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringApply routine.
  SUBROUTINE Petsc_MatColoringApply(matColoring,isColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to apply
    TYPE(PetscISColoringType), INTENT(INOUT) :: isColoring !<The index set coloring to apply to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringApply",err,error,*999)

    CALL MatColoringApply(matColoring%matColoring,isColoring%isColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringApply.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringApply")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringApply",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringApply
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringCreate routine.
  SUBROUTINE Petsc_MatColoringCreate(a,matColoring,err,error,*)

    !Argument Variables
    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to create the coloring for
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<On exit, the created coloring
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringCreate",err,error,*999)

    CALL MatColoringCreate(a%mat,matColoring%matColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringCreate")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringCreate
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringDestroy routine.
  SUBROUTINE Petsc_MatColoringDestroy(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringDestroy",err,error,*999)

    CALL MatColoringDestroy(matColoring%matColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringDestroy.",err,error,*999)
    ENDIF
    matColoring%matColoring=PETSC_NULL_OBJECT
    
    EXITS("Petsc_MatColoringDestroy")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringDestroy",err,error)
    EXITS("Petsc_MatColoringDestroy")
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringDestroy
    
  !
  !================================================================================================================================
  !
  
  !Finalise the PETSc MatColoring structure and destroy the MatColoring
  SUBROUTINE Petsc_MatColoringFinalise(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The matColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringFinalise",err,error,*999)

    IF(matColoring%matColoring/=PETSC_NULL_OBJECT) THEN
      CALL Petsc_MatColoringDestroy(matColoring,err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringFinalise")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc MatColoring structure
  SUBROUTINE Petsc_MatColoringInitialise(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The matColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringInitialise",err,error,*999)

    matColoring%matColoring=PETSC_NULL_OBJECT
    
    EXITS("Petsc_MatColoringInitialise")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringInitialise
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringSetFromOptions routine.
  SUBROUTINE Petsc_MatColoringSetFromOptions(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringSetFromOptions",err,error,*999)

    CALL MatColoringSetFromOptions(matColoring%matColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringSetFromOptions
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringSetType routine.
  SUBROUTINE Petsc_MatColoringSetType(matColoring,coloringType,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to set the type for
    MatColoringType, INTENT(IN) :: coloringType !<The type of coloring to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringSetType",err,error,*999)

    CALL MatColoringSetType(matColoring%matColoring,coloringType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringSetType")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringSetType
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreate routine.
  SUBROUTINE PETSC_MATCREATE(COMMUNICATOR,A,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI Communicator
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATE",err,error,*999)

    CALL MatCreate(COMMUNICATOR,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATE")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATE    
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,3,0) )
  !>Buffer routine to the PETSc MatCreateAIJ routine.
  SUBROUTINE PETSC_MATCREATEAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
    & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A,err,error,*)

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
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATEAIJ",err,error,*999)

    CALL MatCreateAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
      & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateAIJ",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATEAIJ")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATEAIJ",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEAIJ
#else  
  !>Buffer routine to the PETSc MatCreateMPIAIJ routine.
  SUBROUTINE PETSC_MATCREATEMPIAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
    & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A,err,error,*)

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
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATEMPIAIJ",err,error,*999)

    CALL MatCreateMPIAIJ(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,DIAG_NUMBER_NZ_PERROW,DIAG_NUMBER_NZ_EACHROW, &
      & OFFDIAG_NUMBER_NZ_PERROW,OFFDIAG_NUMBER_NZ_EACHROW,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateMPIAIJ",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATEMPIAIJ")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATEMPIAIJ",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEMPIAIJ
#endif  
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,3,0) )
  !>Buffer routine to the PETSc MatCreateDense routine.
  SUBROUTINE PETSC_MATCREATEDENSE(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<The number of global rows
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*) !<Optional, the allocated matrix data.
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATEDENSE",err,error,*999)

    CALL MatCreateDense(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateDense",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATEDENSE")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATEDENSE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEDENSE
#else  
  !>Buffer routine to the PETSc MatCreateMPIDense routine.
  SUBROUTINE PETSC_MATCREATEMPIDENSE(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<The number of global rows
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*) !<Optional, the allocated matrix data.
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the matrix to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATEMPIDENSE",err,error,*999)

    CALL MatCreateMPIDense(COMMUNICATOR,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,MATRIX_DATA,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateMPIDense",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATEMPIDENSE")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATEMPIDENSE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATEMPIDENSE
#endif
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateSeqAIJ routine.
  SUBROUTINE PETSC_MATCREATESEQAIJ(COMMUNICATOR,M,N,NUMBER_NZ_PERROW,NUMBER_NZ_EACHROW,A,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: M !<The number of rows
    INTEGER(INTG), INTENT(IN) :: N !<The number of columns
    INTEGER(INTG), INTENT(IN) :: NUMBER_NZ_PERROW !<The maximum number of non-zeros per row
    INTEGER(INTG), INTENT(IN) :: NUMBER_NZ_EACHROW(*) !<The number of non-zeros in each row
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATESEQAIJ",err,error,*999)

    CALL MatCreateSeqAIJ(COMMUNICATOR,M,N,NUMBER_NZ_PERROW,NUMBER_NZ_EACHROW,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateSeqAIJ",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATESEQAIJ")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATESEQAIJ",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATESEQAIJ
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateSeqDense routine.
  SUBROUTINE PETSC_MATCREATESEQDENSE(COMMUNICATOR,M,N,MATRIX_DATA,A,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI Communicator
    INTEGER(INTG), INTENT(IN) :: M !<The number of rows
    INTEGER(INTG), INTENT(IN) :: N !<The number of columns
    REAL(DP), INTENT(IN) :: MATRIX_DATA(*) !<Optional, the allocated matrix data
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATCREATESEQDENSE",err,error,*999)

    CALL MatCreateSeqDense(COMMUNICATOR,M,N,MATRIX_DATA,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateSeqDense",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATCREATESEQDENSE")
    RETURN
999 ERRORSEXITS("PETSC_MATCREATESEQDENSE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATCREATESEQDENSE

  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatSetType routine.
  SUBROUTINE PETSC_MATSETTYPE(A,MATRIXTYPE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix
    MatType, INTENT(IN) :: MATRIXTYPE !<The matrix type 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETTYPE",err,error,*999)

    CALL MatSetType(A%mat,MATRIXTYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_MATSETTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETTYPE

  !
  !================================================================================================================================
  !    

  !>Buffer routine to the PETSc MatDestroy routine.
  SUBROUTINE PETSC_MATDESTROY(A,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATDESTROY",err,error,*999)

    CALL MatDestroy(A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatDestroy",err,error,*999)
    ENDIF
    A%mat=PETSC_NULL_OBJECT
     
    EXITS("PETSC_MATDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_MATDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_MATDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringCreate routine.
  SUBROUTINE Petsc_MatFDColoringCreate(A,isColoring,matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The PETSc matrix to create the FD coloring for
    TYPE(PetscISColoringType), INTENT(IN) :: isColoring !<The index set coloring to create the finite difference coloring for
    TYPE(PetscMatFDColoringType), INTENT(OUT) :: matFDColoring !<On exit, the matrix finite difference coloring
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringCreate",err,error,*999)

    CALL MatFDColoringCreate(A%mat,isColoring%isColoring,matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringCreate")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringDestroy routine.
  SUBROUTINE Petsc_MatFDColoringDestroy(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringDestroy",err,error,*999)

    CALL MatFDColoringDestroy(matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringDestroy.",err,error,*999)
    ENDIF
    matFDColoring%matFDColoring=PETSC_NULL_OBJECT
     
    EXITS("Petsc_MatFDColoringDestroy")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringDestroy
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc MatFDColoring structure and destroy the MatFDColoring
  SUBROUTINE Petsc_MatFDColoringFinalise(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The MatFDColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringFinalise",err,error,*999)

    IF(matFDColoring%matFDColoring/=PETSC_NULL_OBJECT) THEN
      CALL Petsc_MatFDColoringDestroy(matFDColoring,err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringFinalise")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc MatFDColoring structure
  SUBROUTINE Petsc_MatFDColoringInitialise(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The MatFDColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringInitialise",err,error,*999)

    matFDColoring%matFDColoring=PETSC_NULL_OBJECT
    
    EXITS("Petsc_MatFDColoringInitialise")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringInitialise
    
  !
  !================================================================================================================================2
  !

  !>Buffer routine to the PETSc MatFDColoringSetFromOptions routine.
  SUBROUTINE Petsc_MatFDColoringSetFromOptions(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetFromOptions",err,error,*999)

    CALL MatFDColoringSetFromOptions(matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetFromOptions
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetParameters routine.
  SUBROUTINE Petsc_MatFDColoringSetParameters(matFDColoring,rError,uMin,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    REAL(DP) :: rError !<The square root of the relative error
    REAL(DP) :: uMin !<MatFDColoring umin option
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetParameters",err,error,*999)

    CALL MatFDColoringSetParameters(matFDColoring%matFDColoring,rError,uMin,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetParameters.",err,error,*999)
    ENDIF

    EXITS("Petsc_MatFDColoringSetParameters")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetParameters",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetParameters

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetFunction routine.
  SUBROUTINE Petsc_MatFDColoringSetFunction(matFDColoring,fFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    EXTERNAL fFunction !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetFunction",err,error,*999)

    CALL MatFDColoringSetFunction(matFDColoring%matFDColoring,fFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetFunction.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetFunction")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetFunction",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetFunction
        
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_LT(3,0,0) )
  !>Buffer routine to the PETSc MatFDColoringSetFunctionSNES routine.
  SUBROUTINE Petsc_MatFDColoringSetFunctionSnes(matFDColoring,fFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    EXTERNAL fFunction !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetFunctionSnes",err,error,*999)

    CALL MatFDColoringSetFunctionSNES(matFDColoring%matFDColoring,fFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetFunctionSNES.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetFunctionSnes")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetFunctionSnes",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetFunctionSnes
#endif
        
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetup routine.
  SUBROUTINE Petsc_MatFDColoringSetup(A,isColoring,matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The PETSc matrix to setup the FD coloring for
    TYPE(PetscISColoringType), INTENT(INOUT) :: isColoring !<The index set coloring to setup the finite difference coloring for
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to setup
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetup",err,error,*999)

    CALL MatFDColoringSetup(A%mat,isColoring%isColoring,matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetup.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetup")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetup
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatDenseGetArrayF90 routine.
  SUBROUTINE Petsc_MatDenseGetArrayF90(A,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT), TARGET :: A !<The matrix to get the array for
    REAL(DP), POINTER :: array(:,:) !<On exit, a pointer to the matrix array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatDenseGetArrayF90",err,error,*999)

    IF(ASSOCIATED(array)) THEN
      CALL FlagError("Array is already associated.",err,error,*999)
    ELSE
      CALL MatDenseGetArrayF90(A%mat,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in MatDenseGetArrayF90.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_MatDenseGetArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatDenseGetArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatDenseGetArrayF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSeqAIJGetArrayF90 routine.
  SUBROUTINE Petsc_MatSeqAIJGetArrayF90(a,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT), TARGET :: A !<The matrix to get the array for
    REAL(DP), POINTER :: array(:,:) !<On exit, a pointer to the matrix array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSeqAIJGetArrayF90",err,error,*999)

    IF(ASSOCIATED(ARRAY)) THEN
      CALL FlagError("Array is already associated.",err,error,*999)
    ELSE
      CALL MatSeqAIJGetArrayF90(A%mat,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in MatSeqAIJGetArrayF90.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_MatSeqAIJGetArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatSeqAIJGetArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSeqAIJGetArrayF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSeqAIJGetMaxRowNonzeros routine.
  SUBROUTINE Petsc_MatSeqAIJGetMaxRowNonzeros(A,maxNumberNonZeros,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to get the maximum number of non zeros for
    INTEGER(INTG), INTENT(OUT) :: maxNumberNonZeros!<On exit, the maximum number of non zeros in any row of the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSeqAIJGetMaxRowNonzeros",err,error,*999)

    !CALL MatSeqAIJGetMaxRowNonzeros(A%mat,maxNumberNonZeros,err)
    CALL FlagError("Not implemented.",err,error,*999)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSeqAIJGetMaxRowNonzeros.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSeqAIJGetMaxRowNonzeros")
    RETURN
999 ERRORSEXITS("Petsc_MatSeqAIJGetMaxRowNonzeros",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSeqAIJGetMaxRowNonzeros
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetInfo routine.
  SUBROUTINE PETSC_MATGETINFO(A,FLAG,INFO,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to get the information of
    MatInfoType, INTENT(IN) :: FLAG !<The matrix information collective to return
    MatInfo, INTENT(OUT) :: INFO(MAT_INFO_SIZE) !<On return, the matrix information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATGETINFO",err,error,*999)

    CALL MatGetInfo(A%MAT,FLAG,INFO,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetInfo",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATGETINFO")
    RETURN
999 ERRORSEXITS("PETSC_MATGETINFO",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATGETINFO
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetOwnershipRange routine.
  SUBROUTINE PETSC_MATGETOWNERSHIPRANGE(A,FIRST_ROW,LAST_ROW,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to get the ownership range of
    INTEGER(INTG), INTENT(OUT) :: FIRST_ROW !<On exit, the first row for the matrix
    INTEGER(INTG), INTENT(OUT) :: LAST_ROW !<On exit, the last row for the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATGETOWNERSHIPRANGE",err,error,*999)

    CALL MatGetOwnershipRange(A%MAT,FIRST_ROW,LAST_ROW,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetOwnershipRange",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATGETOWNERSHIPRANGE")
    RETURN
999 ERRORSEXITS("PETSC_MATGETOWNERSHIPRANGE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_MATGETOWNERSHIPRANGE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetRow routine.
  SUBROUTINE Petsc_MatGetRow(A,rowNumber,numberOfColumns,columns,values,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to get the array for
    INTEGER(INTG), INTENT(IN) :: rowNumber !<The row number to get the row values for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On return, the number of nonzero columns in the row
    INTEGER(INTG), INTENT(OUT) :: columns(:) !<On return, the column numbers for the nonzero columns in the row
    REAL(DP), INTENT(OUT) :: values(:) !<On exit, the nonzero values in the row.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatGetRow",err,error,*999)

    CALL MatGetRow(A%mat,rowNumber,numberOfColumns,columns,values,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetRow.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatGetRow")
    RETURN
999 ERRORSEXITS("Petsc_MatGetRow",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatGetRow
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetValues routine.
  SUBROUTINE PETSC_MATGETVALUES(A,M,M_INDICES,N,N_INDICES,VALUES,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to get the values of
    INTEGER(INTG), INTENT(IN) :: M !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: N !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*) !<The column indices
    REAL(DP), INTENT(OUT) :: VALUES(*) !<The values to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATGETVALUES",err,error,*999)

    CALL MatGetValues(A%mat,M,M_INDICES,N,N_INDICES,VALUES,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetValues",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATGETVALUES")
    RETURN
999 ERRORSEXITS("PETSC_MATGETVALUES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATGETVALUES
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc Mat structure and destroy the KSP
  SUBROUTINE PETSC_MATFINALISE(MAT_,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: MAT_ !<The MAT to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATFINALISE",err,error,*999)

    IF(MAT_%mat/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_MATDESTROY(MAT_,err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_MATFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_MATFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Mat structure
  SUBROUTINE PETSC_MATINITIALISE(A,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The MAT to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATINITIALISE",err,error,*999)

    A%mat=PETSC_NULL_OBJECT
    !MAT_%mat_DATA(1)=0
    !MAT_%mat_OFFSET=0
    
    EXITS("PETSC_MATINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_MATINITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_MATINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatDenseRestoreArrayF90 routine.
  SUBROUTINE Petsc_MatDenseRestoreArrayF90(A,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to restore the array for
    REAL(DP), POINTER :: array(:,:) !<A pointer to the matrix array to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatDenseRestoreArrayF90",err,error,*999)

    CALL MatDenseRestoreArrayF90(A%MAT,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatDenseRestoreArrayF90.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatDenseRestoreArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatDenseRestoreArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatDenseRestoreArrayF90
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatSeqAIJRestoreArrayF90 routine.
  SUBROUTINE Petsc_MatSeqAIJRestoreArrayF90(A,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to restore the array for
    REAL(DP), POINTER :: array(:,:) !<A pointer to the matrix array to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MatSeqAIJRestoreArrayF90",err,error,*999)

    CALL MatSeqAIJRestoreArrayF90(A%mat,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSeqAIJRestoreArrayF90.",err,error,*999)
    ENDIF
    
    EXITS("MatSeqAIJRestoreArrayF90")
    RETURN
999 ERRORSEXITS("MatSeqAIJRestoreArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE 
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatRestoreRow routine.
  SUBROUTINE Petsc_MatRestoreRow(A,rowNumber,numberOfColumns,columns,values,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to restore the row for
    INTEGER(INTG), INTENT(IN) :: rowNumber !<The row number to restore the row values for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<The number of nonzero columns in the row
    INTEGER(INTG), INTENT(OUT) :: columns(:) !<The column numbers for the nonzero columns in the row
    REAL(DP), INTENT(OUT) :: values(:) !<The nonzero values in the row to restore.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatRestoreRow",err,error,*999)

    CALL MatRestoreRow(A%mat,rowNumber,numberOfColumns,columns,values,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatRestoreRow.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatRestoreRow")
    RETURN
999 ERRORSEXITS("Petsc_MatRestoreRow",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatRestoreRow
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetLocalToGlobalMapping routine.
  SUBROUTINE PETSC_MATSETLOCALTOGLOBALMAPPING(A,CTX,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the local to global mapping for
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETLOCALTOGLOBALMAPPING",err,error,*999)

    CALL MatSetLocalToGlobalMapping(A%mat,CTX%ISLOCALTOGLOBALMAPPING,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetLocalToGlobalMapping",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETLOCALTOGLOBALMAPPING")
    RETURN
999 ERRORSEXITS("PETSC_MATSETLOCALTOGLOBALMAPPING",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETLOCALTOGLOBALMAPPING
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetOption routine.
  SUBROUTINE PETSC_MATSETOPTION(A,OPTION,FLAG,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the option for
    MatOption, INTENT(IN) :: OPTION !<The option to set
#if ( PETSC_VERSION_GE(3,2,0) )
    PetscBool, INTENT(IN) :: FLAG !<The option flag to set
#else
    PetscTruth, INTENT(IN) :: FLAG !<The option flag to set
#endif    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETOPTION",err,error,*999)

    CALL MatSetOption(A%mat,OPTION,FLAG,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetOption",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETOPTION")
    RETURN
999 ERRORSEXITS("PETSC_MATSETOPTION",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETOPTION
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetSizes routine.
  SUBROUTINE PETSC_MATSETSIZES(A,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the size of
    INTEGER(INTG), INTENT(IN) :: LOCAL_M !<Number of local rows
    INTEGER(INTG), INTENT(IN) :: LOCAL_N !<Number of local columns
    INTEGER(INTG), INTENT(IN) :: GLOBAL_M !<Number of global rows
    INTEGER(INTG), INTENT(IN) :: GLOBAL_N !<Number of global columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETSIZES",err,error,*999)

    CALL MatSetSizes(A%mat,LOCAL_M,LOCAL_N,GLOBAL_M,GLOBAL_N,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetSizes",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETSIZES")
    RETURN
999 ERRORSEXITS("PETSC_MATSETSIZES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETSIZES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValue routine.
  SUBROUTINE PETSC_MATSETVALUE(A,ROW,COL,VALUE,INSERT_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: ROW !<The row index
    INTEGER(INTG), INTENT(IN) :: COL !<The column index
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETVALUE",err,error,*999)

    CALL MatSetValue(A%mat,ROW,COL,VALUE,INSERT_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValue",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETVALUE")
    RETURN
999 ERRORSEXITS("PETSC_MATSETVALUE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValues routine.
  SUBROUTINE PETSC_MATSETVALUES(A,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: M !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: N !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*) !<The column indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETVALUES",err,error,*999)

    CALL MatSetValues(A%mat,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValues",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETVALUES")
    RETURN
999 ERRORSEXITS("PETSC_MATSETVALUES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValueLocal routine.
  SUBROUTINE PETSC_MATSETVALUELOCAL(A,ROW,COL,VALUE,INSERT_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: ROW !<The row index
    INTEGER(INTG), INTENT(IN) :: COL !<The column index
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETVALUELOCAL",err,error,*999)

    CALL MatSetValueLocal(A%mat,ROW,COL,VALUE,INSERT_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValueLocal",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETVALUELOCAL")
    RETURN
999 ERRORSEXITS("PETSC_MATSETVALUELOCAL",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUELOCAL
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValuesLocal routine.
  SUBROUTINE PETSC_MATSETVALUESLOCAL(A,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: M !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: M_INDICES(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: N !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: N_INDICES(*) !<The column indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATSETVALUESLOCAL",err,error,*999)

    CALL MatSetValuesLocal(A%mat,M,M_INDICES,N,N_INDICES,VALUES,INSERT_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValuesLocal",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATSETVALUESLOCAL")
    RETURN
999 ERRORSEXITS("PETSC_MATSETVALUESLOCAL",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATSETVALUESLOCAL
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatView routine.
  SUBROUTINE PETSC_MATVIEW(A,V,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to view
    PetscViewer, INTENT(IN) :: V !<The viewer
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATVIEW",err,error,*999)

    CALL MatView(A%mat,V,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatView",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATVIEW")
    RETURN
999 ERRORSEXITS("PETSC_MATVIEW",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATVIEW
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatZeroEntries routine.
  SUBROUTINE PETSC_MATZEROENTRIES(A,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to zero the entries of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_MATZEROENTRIES",err,error,*999)

    CALL MatZeroEntries(A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatZeroEntries",err,error,*999)
    ENDIF
    
    EXITS("PETSC_MATZEROENTRIES")
    RETURN
999 ERRORSEXITS("PETSC_MATZEROENTRIES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_MATZEROENTRIES
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc PC structure and destroy the KSP
  SUBROUTINE PETSC_PCFINALISE(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The PC to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_PCFINALISE",err,error,*999)

    IF(pc%pc/=PETSC_NULL_OBJECT) THEN
      !Do nothing - should be destroyed when the KSP is destroyed.
    ENDIF
    
    EXITS("PETSC_PCFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_PCFINALISE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_PCFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc PC structure
  SUBROUTINE PETSC_PCINITIALISE(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The PC to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_PCINITIALISE",err,error,*999)

    pc%pc=PETSC_NULL_OBJECT
    
    EXITS("PETSC_PCINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_PCINITIALISE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_PCINITIALISE
    
  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,0,0) )
  !>Buffer routine to the PETSc PCFactoSetMatSolverPackage routine.
  SUBROUTINE PETSC_PCFACTORSETMATSOLVERPACKAGE(pc,SOLVER_PACKAGE,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the solver package for
    MatSolverPackage, INTENT(IN) :: SOLVER_PACKAGE !<The solver package to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_PCFACTORSETMATSOLVERPACKAGE",err,error,*999)

    CALL PCFactorSetMatSolverPackage(pc%pc,SOLVER_PACKAGE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCFactorSetMatSolverPackage",err,error,*999)
    ENDIF
    
    EXITS("PETSC_PCFACTORSETMATSOLVERPACKAGE")
    RETURN
999 ERRORSEXITS("PETSC_PCFACTORSETMATSOLVERPACKAGE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_PCFACTORSETMATSOLVERPACKAGE
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,0,0) )
  !>Buffer routine to the PETSc PCFactorSetUpMatSolverPackage routine.
  SUBROUTINE Petsc_PCFactorSetUpMatSolverPackage(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the solver package for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCFactorSetUpMatSolverPackage",err,error,*999)

    CALL PCFactorSetUpMatSolverPackage(pc%pc,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCFactorSetUpMatSolverPackage",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCFactorSetUpMatSolverPackage")
    RETURN
999 ERRORSEXITS("Petsc_PCFactorSetUpMatSolverPackage",err,error)
    RETURN 1
  END SUBROUTINE Petsc_PCFactorSetUpMatSolverPackage
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,0,0) )
  !>Buffer routine to the PETSc PCFactorGetMatrix routine.
  SUBROUTINE Petsc_PCFactorGetMatrix(pc,factoredMatrix,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the solver package for
    TYPE(PetscMatType), INTENT(OUT) :: factoredMatrix !<The factored matrix to get from preconditioner context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCFactorGetMatrix",err,error,*999)

    CALL PCFactorGetMatrix(pc%pc,factoredMatrix%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCFactorGetMatrix",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCFactorGetMatrix")
    RETURN
999 ERRORSEXITS("Petsc_PCFactorGetMatrix",err,error)
    RETURN 1
  END SUBROUTINE Petsc_PCFactorGetMatrix
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,0,0) )
  !>Buffer routine to the PETSc MatMumpsSetIcntl routine.
  SUBROUTINE Petsc_MatMumpsSetIcntl(factoredMatrix,icntl,ival,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: factoredMatrix !<The factored matrix from PETSc-MUMPS interface
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS ICNTL integer control parameter
    INTEGER(INTG), INTENT(IN) :: ival !<The MUMPS ICNTL integer value to set: ICNTL(icntl)=ival
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatMumpsSetIcntl",err,error,*999)

    CALL MatMumpsSetIcntl(factoredMatrix%mat,icntl,ival,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatMumpsSetIcntl",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatMumpsSetIcntl")
    RETURN
999 ERRORSEXITS("Petsc_MatMumpsSetIcntl",err,error)
    RETURN 1
  END SUBROUTINE Petsc_MatMumpsSetIcntl
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,4,0) )
  !>Buffer routine to the PETSc MatMumpsSetCntl routine.
  SUBROUTINE Petsc_MatMumpsSetCntl(factoredMatrix,icntl,val,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: factoredMatrix !<The factored matrix from PETSc-MUMPS interface
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS CNTL integer control parameter
    REAL(DP), INTENT(IN) :: val !<The MUMPS CNTL real value to set: CNTL(icntl)=val
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatMumpsSetCntl",err,error,*999)

    CALL MatMumpsSetCntl(factoredMatrix%mat,icntl,val,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatMumpsSetCntl.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatMumpsSetCntl")
    RETURN
999 ERRORSEXITS("Petsc_MatMumpsSetCntl",err,error)
    RETURN 1
  END SUBROUTINE Petsc_MatMumpsSetCntl
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )
  !>Buffer routine to the PETSc PCSetReusePreconditioner routine.
  SUBROUTINE PETSC_PCSETREUSEPRECONDITIONER(pc,FLAG,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the type of
    PetscBool, INTENT(IN) :: FLAG !<True/false
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_PCSETREUSEPRECONDITIONER",err,error,*999)

    CALL PCSetReusePreconditioner(pc%pc,FLAG,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCSetReusePreconditioner",err,error,*999)
    ENDIF
    
    EXITS("PETSC_PCSETREUSEPRECONDITIONER")
    RETURN
999 ERRORSEXITS("PETSC_PCSETREUSEPRECONDITIONER",err,error)
    RETURN 1
  END SUBROUTINE PETSC_PCSETREUSEPRECONDITIONER
#endif
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCSetType routine.
  SUBROUTINE PETSC_PCSETTYPE(pc,METHOD,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the type of
    PCType, INTENT(IN) :: METHOD !<The preconditioning method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_PCSETTYPE",err,error,*999)

    CALL PCSetType(pc%pc,METHOD,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCSetType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_PCSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_PCSETTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_PCSETTYPE
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc SNES structure and destroy the SNES
  SUBROUTINE PETSC_SNESFINALISE(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESFINALISE",err,error,*999)

    IF(snes%snes/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_SNESDESTROY(snes,err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_SNESFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_SNESFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc SNES structure
  SUBROUTINE PETSC_SNESINITIALISE(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The snes to 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESINITIALISE",err,error,*999)

    snes%snes=PETSC_NULL_OBJECT
     
    EXITS("PETSC_SNESINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_SNESINITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_SNESINITIALISE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESCreate routine.
  SUBROUTINE PETSC_SNESCREATE(COMMUNICATOR,snes,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator for the SNES creation
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<On exit, the SNES information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESCREATE",err,error,*999)

    CALL SNESCreate(COMMUNICATOR,snes%snes,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESCreate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESCREATE")
    RETURN
999 ERRORSEXITS("PETSC_SNESCREATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESDestroy routine.
  SUBROUTINE PETSC_SNESDESTROY(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESDESTROY",err,error,*999)

    CALL SNESDestroy(snes%snes,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESDestroy",err,error,*999)
    ENDIF
    snes%snes=PETSC_NULL_OBJECT
     
    EXITS("PETSC_SNESDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_SNESDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_SNESDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetApplicationContext routine.
  SUBROUTINE Petsc_SnesGetApplicationContext(snes,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the context for
    TYPE(SOLVER_TYPE), POINTER :: ctx !<On exit, the solver data context to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetApplicationContext",err,error,*999)

    IF(ASSOCIATED(ctx)) THEN
      CALL FlagError("Context is already associated.",err,error,*999)
    ELSE
      CALL SNESGetApplicationContext(snes%snes,ctx,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in SNESGetApplicationContext.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_SnesGetApplicationContext")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetApplicationContext",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetApplicationContext

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetConvergedReason routine.
  SUBROUTINE PETSC_SNESGETCONVERGEDREASON(snes,REASON,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the converged reason for
    INTEGER(INTG), INTENT(OUT) :: REASON !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESGETCONVERGEDREASON",err,error,*999)

    CALL SNESGetConvergedReason(snes%snes,REASON,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetConvergedReason",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESGETCONVERGEDREASON")
    RETURN
999 ERRORSEXITS("PETSC_SNESGETCONVERGEDREASON",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESGETCONVERGEDREASON
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSC SNESGetSolutionUpdate routine.
  SUBROUTINE Petsc_SnesGetSolutionUpdate(snes_,solutionUpdate,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes_ !<The SNES to get the solution update for
    TYPE(PetscVecType), INTENT(INOUT) :: solutionUpdate !<On exit, the solution update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetSolutionUpdate",err,error,*999)

    CALL SNESGetSolutionUpdate(snes_%snes,solutionUpdate%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetSolutionUpdate",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesGetSolutionUpdate")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetSolutionUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetSolutionUpdate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetApplicationContext routine.
  SUBROUTINE Petsc_SnesSetApplicationContext(snes,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the context for
    TYPE(SOLVER_TYPE), POINTER :: ctx !<The solver data context to set as a context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetApplicationContext",err,error,*999)

    CALL SNESSetApplicationContext(snes%snes,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetApplicationContext",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetApplicationContext")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetApplicationContext",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetApplicationContext

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_LT(3,6,0) )
  !>Buffer routine to the PETSc SNESSetFunctionNorm routine.
  SUBROUTINE PETSC_SNESSETFUNCTIONNORM(snes,FUNCTION_NORM,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the function norm for
    REAL(DP), INTENT(OUT) :: FUNCTION_NORM !<On exit, the function norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETFUNCTIONNORM",err,error,*999)

    CALL SNESSetFunctionNorm(snes%snes,FUNCTION_NORM,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetFunctionNorm",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETFUNCTIONNORM")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETFUNCTIONNORM",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETFUNCTIONNORM
#endif    

  !
  !================================================================================================================================
  !

  !>Buffer routine to the petsc SnesLineSearchSetNorms routine.
  SUBROUTINE PETSC_SnesLineSearchSetNorms(snes,XNORM,FNORM,YNORM,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the computed norms for X, Y, and F
    REAL(DP), INTENT(INOUT) :: XNORM !<On exit, the norm of the current solution
    REAL(DP), INTENT(INOUT) :: FNORM !<On exit, the norm of the current function
    REAL(DP), INTENT(INOUT) :: YNORM !<On exit, the norm of the current update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("petsc_SnesLineSearchSetNorms",err,error,*999)

    CALL SnesLineSearchSetNorms(snes%snes,XNORM,FNORM,YNORM,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("petsc error in SnesLineSearchSetNorms",err,error,*999)
    ENDIF
    
    EXITS("petsc_SnesLineSearchSetNorms")
    RETURN
999 ERRORSEXITS("petsc_SnesLineSearchSetNorms",err,error)
    RETURN 1
  END SUBROUTINE petsc_SnesLineSearchSetNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the petsc SnesLineSearchGetNorms routine.
  SUBROUTINE petsc_SnesLineSearchGetNorms(lineSearch,XNORM,FNORM,YNORM,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to get the norms for X, Y, and F from.
    REAL(DP), INTENT(INOUT) :: XNORM !<On exit, the norm of the current solution
    REAL(DP), INTENT(INOUT) :: FNORM !<On exit, the norm of the current function
    REAL(DP), INTENT(INOUT) :: YNORM !<On exit, the norm of the current update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("petsc_SnesLineSearchGetNorms",err,error,*999)

    CALL SnesLineSearchGetNorms(lineSearch%snesLineSearch,XNORM,FNORM,YNORM,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("petsc error in SnesLineSearchGetNorms",err,error,*999)
    ENDIF
    
    EXITS("petsc_SnesLineSearchGetNorms")
    RETURN
999 ERRORSEXITS("petsc_SnesLineSearchGetNorms",err,error)
    RETURN 1
  END SUBROUTINE petsc_SnesLineSearchGetNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetIterationNumber routine.
  SUBROUTINE PETSC_SNESGETITERATIONNUMBER(snes,ITERATION_NUMBER,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the iteration number for
    INTEGER(INTG), INTENT(OUT) :: ITERATION_NUMBER !<On exit, the number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESGETITERATIONNUMBER",err,error,*999)

    CALL SNESGetIterationNumber(snes%snes,ITERATION_NUMBER,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetIterationNumber",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESGETITERATIONNUMBER")
    RETURN
999 ERRORSEXITS("PETSC_SNESGETITERATIONNUMBER",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESGETITERATIONNUMBER
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetKSP routine.
  SUBROUTINE PETSC_SNESGETKSP(snes,ksp,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the iteration number for
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<On exit, the KSP to associated with the SNES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESGETKSP",err,error,*999)

    CALL SNESGetKSP(snes%snes,ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetKSP",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESGETKSP")
    RETURN
999 ERRORSEXITS("PETSC_SNESGETKSP",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESGETKSP

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetSNESLineSearch routine.
  SUBROUTINE Petsc_SnesGetSnesLineSearch(snes,lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the SNES line search for
    TYPE(PetscSnesLinesearchType), INTENT(OUT) :: lineSearch !<On return, the SNES line search object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesGetSnesLineSearch",err,error,*999)
    
#if ( PETSC_VERSION_GE(3,4,0) )
    CALL SNESGetLineSearch(snes%snes,lineSearch%snesLineSearch,err)
#else
    CALL SNESGetSNESLineSearch(snes%snes,lineSearch%snesLineSearch,err)
#endif
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetLineSearch",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesGetSnesLineSearch")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetSnesLineSearch",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetSnesLineSearch
  

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetLineSearch routine.
  SUBROUTINE Petsc_SnesGetLineSearch(snes,lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the SNES line search for
    TYPE(PetscSnesLinesearchType), INTENT(OUT) :: lineSearch !<On return, the SNES line search object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesGetLineSearch",err,error,*999)

    CALL SNESGetLineSearch(snes%snes,lineSearch%snesLineSearch,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetLineSearch",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesGetLineSearch")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetLineSearch",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetLineSearch
  


  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,3,0) )
  !Finalise the PETSc SNES line search structure
  SUBROUTINE Petsc_SnesLineSearchFinalise(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesLineSearchFinalise",err,error,*999)

    ! We don't actually call PETSc's SNESLineSearchDestroy as PETSc accesses
    ! the LineSearch when calling SNESDestroy and also destroys it when
    ! calling SNESDestroy, so we'll just let PETSc clean it up.
   lineSearch%snesLineSearch=PETSC_NULL_OBJECT

    EXITS("Petsc_SnesLineSearchFinalise")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchFinalise",err,error)
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

    ENTERS("Petsc_SnesLineSearchInitialise",err,error,*999)

    lineSearch%snesLineSearch=PETSC_NULL_OBJECT

    EXITS("Petsc_SnesLineSearchInitialise")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchInitialise",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchInitialise
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_LT(3,3,0) )
  !>Buffer routine to the PETSc SNESLineSearchSet routine.
  SUBROUTINE PETSC_SNESLINESEARCHSET(snes,LINESEARCH_TYPE,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the line search for
    INTEGER(INTG), INTENT(IN) :: LINESEARCH_TYPE !<The line search type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESLINESEARCHSET",err,error,*999)

    SELECT CASE(LINESEARCH_TYPE)
    CASE(PETSC_SNES_LINESEARCH_NONORMS)
      CALL SNESLineSearchSet(snes%snes,SNESLINESEARCHNONORMS,PETSC_NULL_OBJECT,err)
    CASE(PETSC_SNES_LINESEARCH_NO)
      CALL SNESLineSearchSet(snes%snes,SNESLINESEARCHNO,PETSC_NULL_OBJECT,err)
    CASE(PETSC_SNES_LINESEARCH_QUADRATIC)
      CALL SNESLineSearchSet(snes%snes,SNESLINESEARCHQUADRATIC,PETSC_NULL_OBJECT,err)
    CASE(PETSC_SNES_LINESEARCH_CUBIC)
      CALL SNESLineSearchSet(snes%snes,SNESLINESEARCHCUBIC,PETSC_NULL_OBJECT,err)      
    CASE DEFAULT
      CALL FlagError("Invalid line search type.",err,error,*999)
    END SELECT
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSet.",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESLINESEARCHSET")
    RETURN
999 ERRORSEXITS("PETSC_SNESLINESEARCHSET",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESLINESEARCHSET
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,2,0) )
  !>Buffer routine to the PETSc SNESLineSearchSetMonitor routine.
  SUBROUTINE Petsc_SnesLineSearchSetMonitor(lineSearch,monitorLinesearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set whether to output linesearch debug information
    PetscBool, INTENT(IN) :: monitorLinesearch !<Whether to output linesearch debug information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchSetMonitor",err,error,*999)

    CALL SnesLineSearchSetMonitor(lineSearch%snesLineSearch,monitorLinesearch,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetMonitor",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetMonitor")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetMonitor",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetMonitor
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,3,0) )
  !>Buffer routine to the PETSc SNESLineSearchSetComputeNorms routine.
  SUBROUTINE Petsc_SnesLineSearchSetComputeNorms(lineSearch,computeNorms,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set whether to compute norms
    PetscBool, INTENT(IN) :: computeNorms !<Whether to compute norms
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchSetComputeNorms",err,error,*999)

    CALL SNESLineSearchSetComputeNorms(lineSearch%snesLineSearch,computeNorms,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetComputeNorms",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetComputeNorms")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetComputeNorms",err,error)
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

    ENTERS("Petsc_SnesLineSearchComputeNorms",err,error,*999)

    CALL SnesLineSearchComputeNorms(lineSearch%snesLineSearch,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SnesLineSearchComputeNorms",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchComputeNorms")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchComputeNorms",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchComputeNorms
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,3,0) )
  !>Buffer routine to the PETSc SNESLineSearchSetOrder routine.
  SUBROUTINE Petsc_SnesLineSearchSetOrder(lineSearch,lineSearchOrder,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the line search order for
    SNESLineSearchOrder, INTENT(IN) :: lineSearchOrder !<The line search order to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchSetOrder",err,error,*999)

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
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetOrder",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetOrder")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetOrder",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetOrder
#endif

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_LT(3,3,0) )
  !>Buffer routine to the PETSc SNESLineSearchSetParams routine.
  SUBROUTINE PETSC_SNESLINESEARCHSETPARAMS(snes,ALPHA,MAXSTEP,STEPTOL,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the line search parameters for
    REAL(DP), INTENT(IN) :: ALPHA !<The scalar such that 0.5f_{n+1} . f_{n+1} <= .5*f_n . f_n - alpha |f_n . J . f_n| 
    REAL(DP), INTENT(IN) :: MAXSTEP !<The maximum norm of the update vector
    REAL(DP), INTENT(IN) :: STEPTOL !<the minimum norm fraction of the the original step after scaling
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESLINESEARCHSETPARAMS",err,error,*999)

    CALL SNESLineSearchSetParams(snes%snes,ALPHA,MAXSTEP,STEPTOL,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetParams",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESLINESEARCHSETPARAMS")
    RETURN
999 ERRORSEXITS("PETSC_SNESLINESEARCHSETPARAMS",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESLINESEARCHSETPARAMS
#endif  
    
  !
  !================================================================================================================================
  !
  
#if ( PETSC_VERSION_GE(3,3,0) )
  !>Buffer routine to the PETSc SNESLineSearchSetType routine.
  SUBROUTINE Petsc_SnesLineSearchSetType(lineSearch,lineSearchType,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the line search type for
    SNESLineSearchType, INTENT(IN) :: lineSearchType !<The line search type to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesLineSearchSetType",err,error,*999)

    CALL SNESLineSearchSetType(LINESEARCH%snesLineSearch,lineSearchType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetType",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetType")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetType",err,error)
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

    ENTERS("Petsc_SnesLineSearchBTSetAlpha",err,error,*999)

    CALL SNESLineSearchBTSetAlpha(lineSearch%snesLineSearch,alpha,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      END IF
      CALL FlagError("PETSc error in SNESLineSearchBTSetAlpha",err,error,*999)
    END IF

    EXITS("Petsc_SnesLineSearchBTSetAlpha")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchBTSetAlpha",err,error)
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

    ENTERS("Petsc_SnesLineSearchSetTolerances",err,error,*999)

    CALL SNESLineSearchSetTolerances(lineSearch%snesLineSearch, &
      & steptol,maxstep,rtol,atol,ltol,maxIt,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      END IF
      CALL FlagError("PETSc error in SNESLineSearchSetTolerances",err,error,*999)
    END IF

    EXITS("Petsc_SnesLineSearchSetTolerances")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetTolerances",err,error)
    RETURN 1

  END SUBROUTINE Petsc_SnesLineSearchSetTolerances
#endif

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESMonitorSet routine.
  SUBROUTINE PETSC_SNESMONITORSET(snes,MFUNCTION,CTX,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set from the command line options
    EXTERNAL :: MFUNCTION !<The external monitor function to set
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESMONITORSET",err,error,*999)

    CALL SNESMonitorSet(snes%snes,MFUNCTION,CTX,PETSC_NULL_FUNCTION,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESMonitorSet",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESMONITORSET")
    RETURN
999 ERRORSEXITS("PETSC_SNESMONITORSET",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESMONITORSET

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5 )   
  !>Buffer routine to the PETSc SNESQNSetType routine.
  SUBROUTINE PETSC_SNESQNSETTYPE(snes,QTYPE,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESQNType, INTENT(IN) :: QTYPE !<The SNES QN type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESQNSETTYPE",err,error,*999)

    CALL SNESQNSetType(snes%snes,QTYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESQNSetType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESQNSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_SNESQNSETTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESQNSETTYPE
#endif

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESQNSetRestartType routine.
  SUBROUTINE PETSC_SNESQNSETRESTARTTYPE(snes,RTYPE,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESQNRestartType, INTENT(IN) :: RTYPE !<The SNES QN restart type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESQNSETRESTARTTYPE",err,error,*999)

    CALL SNESQNSetRestartType(snes%snes,RTYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESQNSetRestartType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESQNSETRESTARTTYPE")
    RETURN
999 ERRORSEXITS("PETSC_SNESQNSETRESTARTTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESQNSETRESTARTTYPE

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESQNSetScaleType routine.
  SUBROUTINE PETSC_SNESQNSETSCALETYPE(snes,STYPE,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESQNScaleType, INTENT(IN) :: STYPE !<The SNES QN scaling type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESQNSETSCALETYPE",err,error,*999)

    CALL SNESQNSetScaleType(snes%snes,STYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESQNSetScaleType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESQNSETSCALETYPE")
    RETURN
999 ERRORSEXITS("PETSC_SNESQNSETSCALETYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESQNSETSCALETYPE

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFromOptions routine.
  SUBROUTINE PETSC_SNESSETFROMOPTIONS(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set from the command line options
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETFROMOPTIONS",err,error,*999)

    CALL SNESSetFromOptions(snes%snes,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetFromOptions",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETFROMOPTIONS")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETFROMOPTIONS",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetFunction routine.
  SUBROUTINE Petsc_SnesGetFunction(snes,f,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the function from
    TYPE(PetscVecType), INTENT(OUT) :: f !<The residual vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetFunction",err,error,*999)

    CALL SNESGetFunction(snes%snes,f%vec,PETSC_NULL_FUNCTION,PETSC_NULL_OBJECT,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetFunction.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetFunction")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetFunction",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetFunction

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunction routine.
  SUBROUTINE PETSC_SNESSETFUNCTION(snes,F,FFUNCTION,CTX,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscVecType), INTENT(INOUT) :: F !<The residual vector
    EXTERNAL FFUNCTION !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETFUNCTION",err,error,*999)

    CALL SNESSetFunction(snes%snes,F%vec,FFUNCTION,CTX,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetFunction",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETFUNCTION")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETFUNCTION",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETFUNCTION

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunction routine.
  SUBROUTINE PETSC_SNESSETCONVERGENCETEST(snes,CFUNCTION,CTX,err,error,*)
    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    EXTERNAL CFUNCTION !<The external function to call (OpenCMISS subroutine to calculate convergence
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the convergence test function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETCONVERGENCETEST",err,error,*999)

    CALL SNESSetConvergenceTest(snes%snes,CFUNCTION,CTX,PETSC_NULL_FUNCTION,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetConvergenceTest",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETCONVERGENCETEST")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETCONVERGENCETEST",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETCONVERGENCETEST


  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for MatFDColoring contexts.
  SUBROUTINE PETSC_SNESSETJACOBIAN_MATFDCOLORING(snes,A,B,jFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The snes to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The Jacobian matrix
    TYPE(PetscMatType), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
    EXTERNAL jFunction !<The external function to call
    TYPE(PetscMatFDColoringType) :: ctx !<The MatFDColoring data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETJACOBIAN_MATFDCOLORING",err,error,*999)

    CALL SNESSetJacobianBuffer(snes,A,B,jFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetJacobian.",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETJACOBIAN_MATFDCOLORING")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETJACOBIAN_MATFDCOLORING",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_SNESSETJACOBIAN_MATFDCOLORING
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESGETJACOBIAN_SOLVER(snes,A,B,JFUNCTION,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The Jacobian matrix
    TYPE(PetscMatType), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
    EXTERNAL JFUNCTION !<The external function to call
!     TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESGETJACOBIAN_SOLVER",err,error,*999)

    CALL SNESGetJacobian(snes%snes,A%mat,B%mat,JFUNCTION,PETSC_NULL_INTEGER,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetJacobian",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESGETJACOBIAN_SOLVER")
    RETURN
999 ERRORSEXITS("PETSC_SNESGETJACOBIAN_SOLVER",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESGETJACOBIAN_SOLVER

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESGETJACOBIAN_SPECIAL(snes,A,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The Jacobian matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESGETJACOBIAN_SPECIAL",err,error,*999)

    CALL SNESGetJacobian(snes%snes,A%mat,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,err)

    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetJacobian",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESGETJACOBIAN_SPECIAL")
    RETURN
999 ERRORSEXITS("PETSC_SNESGETJACOBIAN_SPECIAL",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESGETJACOBIAN_SPECIAL

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE PETSC_SNESSETJACOBIAN_SOLVER(snes,A,B,JFUNCTION,CTX,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The Jacobian matrix
    TYPE(PetscMatType), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
    EXTERNAL JFUNCTION !<The external function to call
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETJACOBIAN_SOLVER",err,error,*999)

    CALL SNESSetJacobian(snes%snes,A%mat,B%mat,JFUNCTION,CTX,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetJacobian",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETJACOBIAN_SOLVER")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETJACOBIAN_SOLVER",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETJACOBIAN_SOLVER

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESComputeJacobianDefaultColor routine.
  SUBROUTINE Petsc_SnesComputeJacobianDefaultColor(snes,x,J,B,ctx,err,error,*)

    !Argument variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc X Vec
    TYPE(PetscMatType), INTENT(INOUT) :: J !<The PETSc J Mat
    TYPE(PetscMatType), INTENT(INOUT) :: B !<The PETSc B Mat
    TYPE(PetscMatFDColoringType), POINTER :: ctx !<The passed through context
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesComputeJacobianDefaultColor",err,error,*999)
    
    CALL SNESComputeJacobianDefaultColor(snes%snes,x%vec,J%mat,B%mat,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESComputeJacobianDefaultColor.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesComputeJacobianDefaultColor")
    RETURN
999 ERRORSEXITS("Petsc_SnesComputeJacobianDefaultColor",err,error)
    RETURN
    
  END SUBROUTINE Petsc_SnesComputeJacobianDefaultColor

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESComputeJacobianDefault routine.
  SUBROUTINE Petsc_SnesComputeJacobianDefault(snes,x,J,B,ctx,err,error,*)

    !Argument variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec
    TYPE(PetscMatType), INTENT(INOUT) :: J !<The PETSc J Mat
    TYPE(PetscMatType), INTENT(INOUT) :: B !<The PETSc B Mat
    TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesComputeJacobianDefault",err,error,*999)
    
    CALL SNESComputeJacobianDefault(snes%snes,X%vec,J%mat,B%mat,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESComputeJacobianDefault",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesComputeJacobianDefault")    
    RETURN
999 ERRORSEXITS("Petsc_SnesComputeJacobianDefault",err,error)
    RETURN
    
  END SUBROUTINE Petsc_SnesComputeJacobianDefault

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetKSP routine.
  SUBROUTINE PETSC_SNESSETKSP(snes,ksp,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the KSP for
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to be associated with the SNES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETKSP",err,error,*999)

    CALL SNESSetKSP(snes%snes,ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetKSP",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETKSP")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETKSP",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_SNESSETKSP

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetTolerances routine.
  SUBROUTINE PETSC_SNESSETTOLERANCES(snes,ABSTOL,RTOL,STOL,MAXIT,MAXF,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the tolerances for
    REAL(DP), INTENT(IN) :: ABSTOL !<The absolute convergence tolerance
    REAL(DP), INTENT(IN) :: RTOL !<The relative convergence tolerance
    REAL(DP), INTENT(IN) :: STOL !<The convergence tolerance for the change in the solution between steps
    INTEGER(INTG), INTENT(IN) :: MAXIT !<The maximum number of iterations
    INTEGER(INTG), INTENT(IN) :: MAXF !<The maximum number of function evaluations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETTOLERANCES",err,error,*999)

    CALL SNESSetTolerances(snes%snes,ABSTOL,RTOL,STOL,MAXIT,MAXF,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetTolerances",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETTOLERANCES")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETTOLERANCES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETTOLERANCES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetTrustRegionTolerance routine.
  SUBROUTINE PETSC_SNESSETTRUSTREGIONTOLERANCE(snes,TRTOL,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the tolerances for
    REAL(DP), INTENT(IN) :: TRTOL !<The trust region tolerance
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETTRUSTREGIONTOLERANCE",err,error,*999)

    CALL SNESSetTrustRegionTolerance(snes%snes,TRTOL,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetTrustRegionTolerance",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETTRUSTREGIONTOLERANCE")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETTRUSTREGIONTOLERANCE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETTRUSTREGIONTOLERANCE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetType routine.
  SUBROUTINE PETSC_SNESSETTYPE(snes,METHOD,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESType, INTENT(IN) :: METHOD !<The SNES type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETTYPE",err,error,*999)

    CALL SNESSetType(snes%snes,METHOD,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETTYPE

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchGetVecs routine.
  SUBROUTINE Petsc_SnesLineSearchGetVecs(lineSearch,x,f,y,w,g,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The PetcsSnesLineSearch to get the vectors from the SNESLineSearch
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The The old solution 
    TYPE(PetscVecType), INTENT(INOUT) :: f !<The old function 
    TYPE(PetscVecType), INTENT(INOUT) :: y !<The search direction 
    TYPE(PetscVecType), INTENT(INOUT) :: w !<The new solution 
    TYPE(PetscVecType), INTENT(INOUT) :: g !<The new function 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchGetVecs",err,error,*999)

    CALL SNESLineSearchGetVecs(lineSearch%snesLineSearch,x%vec,f%vec,y%vec,w%vec,g%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchGetVecs",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesLineSearchGetVecs")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchGetVecs",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchGetVecs

  !
  !================================================================================================================================
  !

#if ( PETSC_VERSION_GE(3,5,0) )
  !>Buffer routine to the PETSc SNESSetNormSchedule routine.
  SUBROUTINE PETSC_SNESSETNORMSCHEDULE(snes,NORMSCHEDULE,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the norm type for
    SNESNormSchedule, INTENT(IN) :: NORMSCHEDULE !<The norm schedule
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETNORMSCHEDULE",err,error,*999)

    CALL SNESSetNormSchedule(snes%snes,NORMSCHEDULE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetNormSchedule",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETNORMSCHEDULE")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETNORMSCHEDULE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_SNESSETNORMSCHEDULE

  !
  !================================================================================================================================
  !
#else

  !>Buffer routine to the PETSc SNESSetNormType routine.
  SUBROUTINE PETSC_SNESSETNORMTYPE(snes,NORMTYPE,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the norm type for
    INTEGER(INTG), INTENT(IN) :: NORMTYPE !<The norm type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSETNORMTYPE",err,error,*999)

    CALL SNESSetNormType(snes%snes,NORMTYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetNormType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSETNORMTYPE")
    RETURN
999 ERRORSEXITS("PETSC_SNESSETNORMTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSETNORMTYPE
#endif  

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSolve routine.
  SUBROUTINE PETSC_SNESSOLVE(snes,B,X,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to solve
    TYPE(PetscVecType), INTENT(INOUT) :: B !<The constant part of the equation
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_SNESSOLVE",err,error,*999)

    CALL SNESSolve(snes%snes,B%vec,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSolve",err,error,*999)
    ENDIF
    
    EXITS("PETSC_SNESSOLVE")
    RETURN
999 ERRORSEXITS("PETSC_SNESSOLVE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_SNESSOLVE

  !
  !================================================================================================================================
  !

  !Finalise the PETSc TS structure and destroy the TS
  SUBROUTINE PETSC_TSFINALISE(ts,err,error,*)

    !Argument Variables
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSFINALISE",err,error,*999)

    IF(ts%ts/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_TSDESTROY(ts,err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_TSFINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc TS structure
  SUBROUTINE PETSC_TSINITIALISE(ts,err,error,*)

    !Argument Variables
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSINITIALISE",err,error,*999)

    ts%ts=PETSC_NULL_OBJECT
     
    EXITS("PETSC_TSINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_TSINITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSINITIALISE
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSCreate routine.
  SUBROUTINE PETSC_TSCREATE(COMMUNICATOR,ts,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(INOUT) :: COMMUNICATOR !<The MPI communicator for the TS creation
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<On exit, the TS information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSCREATE",err,error,*999)

    CALL TSCreate(COMMUNICATOR,ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSCreate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSCREATE")
    RETURN
999 ERRORSEXITS("PETSC_TSCREATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSCREATE
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSDestroy routine.
  SUBROUTINE PETSC_TSDESTROY(ts,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSDESTROY",err,error,*999)

    CALL TSDestroy(ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSDestroy",err,error,*999)
    ENDIF
    ts%ts=PETSC_NULL_OBJECT
    
    EXITS("PETSC_TSDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_TSDESTROY",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSDESTROY
    
  !
  !================================================================================================================================
  !
    
!  !>Buffer routine to the PETSc TSGetApplicationContext routine.
!  SUBROUTINE PETSC_TSGETAPPLICATIONCONTEXT(ts,USERP,err,error,*)

!    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to get the application context from
!    TYPE(SOLVER_TYPE), POINTER :: USERP !<On exit, a pointer to the user application context
!    INTEGER(INTG), INTENT(OUT) :: err !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
!    !Local Variables

!    ENTERS("PETSC_TSSETAPPLICATIONCONTEXT",err,error,*999)

!    IF(ASSOCIATED(USERP)) THEN
!      CALL FlagError("User application pointer is already associated.",err,error,*999)
!    ELSE
!      NULLIFY(USERP)
!      CALL TSGetApplicationContext(ts%ts,USERP,err)
!      IF(err/=0) THEN
!        IF(petscHandleError) THEN
!          CHKERRQ(err)
!        ENDIF
!        CALL FlagError("PETSc error in TSGetApplicationContext",err,error,*999)
!      ENDIF
!    ENDIF
    
!    EXITS("PETSC_TSGETAPPLICATIONCONTEXT")
!    RETURN
!999 ERRORSEXITS("PETSC_TSGETAPPLICATIONCONTEXT",err,error)
!    RETURN 1
!  END SUBROUTINE PETSC_TSGETAPPLICATIONCONTEXT
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TSMonitorSet routine.
  SUBROUTINE PETSC_TSMONITORSET(ts,MFUNCTION,CTX,err,error,*)

    !Argument Variables
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the monitor for
    EXTERNAL :: MFUNCTION !<The external monitor function to set
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSMONITORSET",err,error,*999)

    CALL TSMonitorSet(ts%ts,MFUNCTION,CTX,PETSC_NULL_FUNCTION,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSMonitorSet",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSMONITORSET")
    RETURN
999 ERRORSEXITS("PETSC_TSMONITORSET",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSMONITORSET
    
  !
  !================================================================================================================================
  !
    
!  !>Buffer routine to the PETSc TSSetApplicationContext routine.
!  SUBROUTINE PETSC_TSSETAPPLICATIONCONTEXT(ts,USERP,err,error,*)

!    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the application context for
!    TYPE(SOLVER_TYPE), POINTER :: USERP !<A pointer to the user application context
!    INTEGER(INTG), INTENT(OUT) :: err !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
!    !Local Variables

!    ENTERS("PETSC_TSSETAPPLICATIONCONTEXT",err,error,*999)

!    CALL TSSetApplicationContext(ts%ts,USERP,err)
!    IF(err/=0) THEN
!      IF(petscHandleError) THEN
!        CHKERRQ(err)
!      ENDIF
!      CALL FlagError("PETSc error in TSSetApplicationContext",err,error,*999)
!    ENDIF
    
!    EXITS("PETSC_TSSETAPPLICATIONCONTEXT")
!    RETURN
!999 ERRORSEXITS("PETSC_TSSETAPPLICATIONCONTEXT",err,error)
!    EXITS("PETSC_TSSETAPPLICATIONCONTEXT")
!    RETURN 1
!  END SUBROUTINE PETSC_TSSETAPPLICATIONCONTEXT
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetDuration routine.
  SUBROUTINE PETSC_TSSETDURATION(ts,MAX_STEPS,MAX_TIME,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set from the options
    INTEGER(INTG), INTENT(IN) :: MAX_STEPS !<The maximum number of steps to use
    REAL(DP), INTENT(IN) :: MAX_TIME !<The maximum time to iteration to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETDURATION",err,error,*999)

    CALL TSSetDuration(ts%ts,MAX_STEPS,MAX_TIME,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetDuration",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETDURATION")
    RETURN
999 ERRORSEXITS("PETSC_TSSETDURATION",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETDURATION
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetFromOptions routine.
  SUBROUTINE PETSC_TSSETFROMOPTIONS(ts,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set from the options
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETFROMOPTIONS",err,error,*999)

    CALL TSSetFromOptions(ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetFromOptions",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETFROMOPTIONS")
    RETURN
999 ERRORSEXITS("PETSC_TSSETFROMOPTIONS",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETFROMOPTIONS

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetExactFinalTime routine.
  SUBROUTINE PETSC_TSSETEXACTFINALTIME(ts,EFTOPT,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the initial time step for
    LOGICAL, INTENT(IN) :: EFTOPT !<The option for exact final time to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETEXACTFINALTIME",err,error,*999)

    CALL TSSetExactFinalTime(ts%ts,EFTOPT,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetExactFinalTime",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETEXACTFINALTIME")
    RETURN
999 ERRORSEXITS("PETSC_TTSSETEXACTFINALTIME",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSSETEXACTFINALTIME    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetInitialTimeStep routine.
  SUBROUTINE PETSC_TSSETINITIALTIMESTEP(ts,INITIAL_TIME,TIME_STEP,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the initial time step for
    REAL(DP), INTENT(IN) :: INITIAL_TIME !<The initial time to set
    REAL(DP), INTENT(IN) :: TIME_STEP !<The time step to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETINITIALTIMESTEP",err,error,*999)

    CALL TSSetInitialTimeStep(ts%ts,INITIAL_TIME,TIME_STEP,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetInitialTimeStep",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETINITIALTIMESTEP")
    RETURN
999 ERRORSEXITS("PETSC_TSSETINITIALTIMESTEP",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETINITIALTIMESTEP

  !
  !================================================================================================================================
  !
    
#if ( PETSC_VERSION_LT(3,2,0) )
  !>Buffer routine to the PETSc TSSetMatrices routine.
  SUBROUTINE PETSC_TSSETMATRICES(ts,ARHS,RHSFUNCTION,ALHS,LHSFUNCTION,FLAG,CTX,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the problem type for
    TYPE(PetscMatType), INTENT(INOUT) :: ARHS !<The RHS matrix
    EXTERNAL RHSFUNCTION !<The external RHS matrix evaluation function to call
    TYPE(PetscMatType), INTENT(INOUT) :: ALHS !<The LHS matrix
    EXTERNAL LHSFUNCTION !<The external LHS matrix evaluation function to call
    INTEGER(INTG), INTENT(IN) :: FLAG !<The matrices structure flag
    TYPE(SOLVER_TYPE), POINTER :: CTX !<The solver data to pass to the matrix evaluations functions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETMATRICES",err,error,*999)

    CALL TSSetMatrices(ts%ts,ARHS%mat,RHSFUNCTION,ALHS%mat,LHSFUNCTION,FLAG,CTX,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetMatrices",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETMATRICES")
    RETURN
999 ERRORSEXITS("PETSC_TSSETMATRICES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETMATRICES
    
  !
  !================================================================================================================================
  !
#endif
  
  !>Buffer routine to the PETSc TSSetProblemType routine.
  SUBROUTINE PETSC_TSSETPROBLEMTYPE(ts,PROB_TYPE,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the problem type for
    INTEGER(INTG), INTENT(IN) :: PROB_TYPE !<The problem type to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETPROBLEMTYPE",err,error,*999)

    CALL TSSetProblemType(ts%ts,PROB_TYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetProblemType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETPROBLEMTYPE")
    RETURN
999 ERRORSEXITS("PETSC_TSSETPROBLEMTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETPROBLEMTYPE
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetRHSFunction routine.
  SUBROUTINE PETSC_TSSETRHSFUNCTION(ts,PETSC_RATES,RHSFUNCTION,CTX,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the problem type for
    TYPE(PetscVecType), INTENT(INOUT) :: PETSC_RATES
    EXTERNAL RHSFUNCTION !<The external RHS function to call
    TYPE(CELLML_PETSC_CONTEXT_TYPE), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETRHSFUNCTION",err,error,*999)

    CALL TSSetRHSFunction(ts%ts,PETSC_RATES%vec,RHSFUNCTION,CTX,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetRHSFunction",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETRHSFUNCTION")
    RETURN
999 ERRORSEXITS("PETSC_TSSETRHSFUNCTION",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETRHSFUNCTION

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetSolution routine.
  SUBROUTINE PETSC_TSSETSOLUTION(ts,INITIALSOLUTION,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time step for
    TYPE(PetscVecType), INTENT(IN) :: INITIALSOLUTION !<The initial solution to be set for the TS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETSOLUTION",err,error,*999)

    CALL TSSetSolution(ts%ts,INITIALSOLUTION%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetSolution",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETSOLUTION")
    RETURN
999 ERRORSEXITS("PETSC_TSSETSOLUTION",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSSETSOLUTION

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSGetSolution routine.
  SUBROUTINE PETSC_TSGETSOLUTION(ts,CURRENTSOLUTION,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time step for
    TYPE(PetscVecType), INTENT(INOUT) :: CURRENTSOLUTION !<The current solution to be set for the TS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETSOLUTION",err,error,*999)

    CALL TSGetSolution(ts%ts,CURRENTSOLUTION%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSGetSolution",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSGETSOLUTION")
    RETURN
999 ERRORSEXITS("PETSC_TSGETSOLUTION",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSGETSOLUTION
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetTimeStep routine.
  SUBROUTINE PETSC_TSSETTIMESTEP(ts,TIME_STEP,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time step for
    REAL(DP), INTENT(IN) :: TIME_STEP !<The time step to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETTIMESTEP",err,error,*999)

    CALL TSSetTimeStep(ts%ts,TIME_STEP,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetTimeStep",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETTIMESTEP")
    RETURN
999 ERRORSEXITS("PETSC_TSSETTIMESTEP",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETTIMESTEP
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetType routine.
  SUBROUTINE PETSC_TSSETTYPE(ts,METHOD,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the type for
    TSType, INTENT(IN) :: METHOD !<The time stepping method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSETTYPE",err,error,*999)

    CALL TSSetType(ts%ts,METHOD,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetType",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_TSSETTYPE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSETTYPE
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSolve routine.
  SUBROUTINE PETSC_TSSOLVE(ts,X,FINALTIME,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to solve
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The solution vector
    REAL(DP), INTENT(OUT) :: FINALTIME
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSOLVE",err,error,*999)

    CALL TSSolve(ts%ts,X%vec,FINALTIME,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSolve",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSOLVE")
    RETURN
999 ERRORSEXITS("PETSC_TSSOLVE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSOLVE
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSStep routine.
  SUBROUTINE PETSC_TSSTEP(ts,STEPS,PTIME,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to step
    INTEGER(INTG), INTENT(IN) :: STEPS !<The number of iterations until termination
    REAL(DP), INTENT(IN) :: PTIME !<The time until termination
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSTEP",err,error,*999)

    CALL TSStep(ts%ts,STEPS,PTIME,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSStep",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSTEP")
    RETURN
999 ERRORSEXITS("PETSC_TSSTEP",err,error)
    RETURN 1
  END SUBROUTINE PETSC_TSSTEP

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSundialsSetType routine.
  SUBROUTINE PETSC_TSSUNDIALSSETTYPE(ts,SUNDIALSTYPE,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to step
    TSSundialsType, INTENT(IN) :: SUNDIALSTYPE
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSUNDIALSSETTYPE",err,error,*999)

    CALL TSSundialsSetType(ts%ts,SUNDIALSTYPE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSUNDIALSSETTYPE",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSUNDIALSSETTYPE")
    RETURN
999 ERRORSEXITS("PETSC_TSSUNDIALSSETTYPE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSSUNDIALSSETTYPE  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSundialsSetType routine.
  SUBROUTINE PETSC_TSSUNDIALSSETTOLERANCE(ts,ABSTOL,RELTOL,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to step
    REAL(DP), INTENT(IN) :: ABSTOL
    REAL(DP), INTENT(IN) :: RELTOL
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_TSSUNDIALSSETTOLERANCE",err,error,*999)

    CALL TSSundialsSetTolerance(ts%ts,ABSTOL,RELTOL,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSUNDIALSSETTOLERANCE",err,error,*999)
    ENDIF
    
    EXITS("PETSC_TSSUNDIALSSETTOLERANCE")
    RETURN
999 ERRORSEXITS("PETSC_TSSUNDIALSSETTYPE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_TSSUNDIALSSETTOLERANCE
  !
  !================================================================================================================================
  !

  !Finalise the PETSc Vec structure and destroy the KSP
  SUBROUTINE PETSC_VECFINALISE(VEC_,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: VEC_ !<The Vec to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECFINALISE",err,error,*999)

    IF(VEC_%vec/=PETSC_NULL_OBJECT) THEN
      CALL PETSC_VECDESTROY(VEC_,err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECFINALISE")
    RETURN
999 ERRORSEXITS("PETSC_VECFINALISE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECFINALISE
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Vec structure
  SUBROUTINE PETSC_VECINITIALISE(VEC_,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: VEC_ !<The Vec to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECINITIALISE",err,error,*999)

    VEC_%vec=PETSC_NULL_OBJECT
    !VEC_%vec_DATA(1)=0
    !VEC_%vec_OFFSET=0
    
    EXITS("PETSC_VECINITIALISE")
    RETURN
999 ERRORSEXITS("PETSC_VECINITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECINITIALISE
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecAssemblyBegin routine.
  SUBROUTINE PETSC_VECASSEMBLYBEGIN(X,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to begin the assembly of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECASSEMBLYBEGIN",err,error,*999)

    CALL VecAssemblyBegin(X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecAssemblyBegin",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECASSEMBLYBEGIN")
    RETURN
999 ERRORSEXITS("PETSC_VECASSEMBLYBEGIN",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECASSEMBLYBEGIN
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecAssemblyEnd routine.
  SUBROUTINE PETSC_VECASSEMBLYEND(X,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to end the assembly of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECASSEMBLYEND",err,error,*999)

    CALL VecAssemblyEnd(X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecAssemblyEnd",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECASSEMBLYEND")
    RETURN
999 ERRORSEXITS("PETSC_VECASSEMBLYEND",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECASSEMBLYEND
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCopy routine.
  SUBROUTINE PETSC_VECCOPY(X,Y,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to copy from
    TYPE(PetscVecType), INTENT(INOUT) :: Y !<The vector to copy to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCOPY",err,error,*999)

    CALL VecCopy(X%vec,Y%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCopy",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCOPY")
    RETURN
999 ERRORSEXITS("PETSC_VECCOPY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECCOPY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreate routine.
  SUBROUTINE PETSC_VECCREATE(COMMUNICATOR,X,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATE",err,error,*999)

    CALL VecCreate(COMMUNICATOR,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATE")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateGhost routine.
  SUBROUTINE PETSC_VECCREATEGHOST(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,X,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    INTEGER(INTG), INTENT(IN) :: NUMBER_GHOST !<The number of ghost elements
    INTEGER(INTG), INTENT(IN) :: GHOSTS(*) !<The global location of the each ghost element
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATEGHOST",err,error,*999)

    CALL VecCreateGhost(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateGhost",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATEGHOST")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATEGHOST",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEGHOST
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateGhostWithArray routine.
  SUBROUTINE PETSC_VECCREATEGHOSTWITHARRAY(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,ARRAY,X,err,error,*)

   !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    INTEGER(INTG), INTENT(IN) :: NUMBER_GHOST !<The number of ghost elements
    INTEGER(INTG), INTENT(IN) :: GHOSTS(*) !<The global location of the each ghost element
    REAL(DP), INTENT(OUT) :: ARRAY(*) !<The preallocated array of matrix data
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATEGHOSTWITHARRAY",err,error,*999)

    CALL VecCreateGhostWithArray(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,NUMBER_GHOST,GHOSTS,ARRAY,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateGhostWithArray",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATEGHOSTWITHARRAY")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATEGHOSTWITHARRAY",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEGHOSTWITHARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateMPI routine.
  SUBROUTINE PETSC_VECCREATEMPI(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,X,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATEMPI",err,error,*999)

    CALL VecCreateMPI(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateMPI",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATEMPI")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATEMPI",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEMPI
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateMPIWithArray routine.
  SUBROUTINE PETSC_VECCREATEMPIWITHARRAY(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,ARRAY,X,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    REAL(DP), INTENT(OUT) :: ARRAY(*) !<The preallocated array for the vector data
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATEMPIWITHARRAY",err,error,*999)

    CALL VecCreateMPIWithArray(COMMUNICATOR,LOCAL_SIZE,GLOBAL_SIZE,ARRAY,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateMPIWithArray",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATEMPIWITHARRAY")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATEMPIWITHARRAY",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATEMPIWITHARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateSeq routine.
  SUBROUTINE PETSC_VECCREATESEQ(COMMUNICATOR,SIZE,X,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: SIZE !<The size of the vector
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATESEQ",err,error,*999)

    CALL VecCreateSeq(COMMUNICATOR,SIZE,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateSeq",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATESEQ")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATESEQ",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATESEQ
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateSeqWithArray routine.
  SUBROUTINE PETSC_VECCREATESEQWITHARRAY(COMMUNICATOR,SIZE,ARRAY,X,err,error,*)

    !Argument Variables
    MPI_Comm, INTENT(IN) :: COMMUNICATOR !<The MPI communicator
    INTEGER(INTG), INTENT(IN) :: SIZE !<The size of the vector
    REAL(DP), INTENT(OUT) :: ARRAY(*) !<The preallocated array for the vector data
    TYPE(PetscVecType), INTENT(INOUT) :: X !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECCREATESEQWITHARRAY",err,error,*999)

    CALL VecCreateSeqWithArray(COMMUNICATOR,SIZE,ARRAY,X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateSeqWithArray",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECCREATESEQWITHARRAY")
    RETURN
999 ERRORSEXITS("PETSC_VECCREATESEQWITHARRAY",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECCREATESEQWITHARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDestroy routine.
  SUBROUTINE PETSC_VECDESTROY(X,err,error,*)

   !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECDESTROY",err,error,*999)

    CALL VecDestroy(X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecDestroy",err,error,*999)
    ENDIF
    X%vec=PETSC_NULL_OBJECT
    
    EXITS("PETSC_VECDESTROY")
    RETURN
999 ERRORSEXITS("PETSC_VECDESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECDESTROY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDuplicate routine.
  SUBROUTINE PETSC_VECDUPLICATE(OLD,NEW,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: OLD !<The vector to duplicate
    TYPE(PetscVecType), INTENT(OUT) :: NEW !<On exit, the new duplicated vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECDUPLICATE",err,error,*999)

    CALL VecDuplicate(OLD%vec,NEW%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecDuplicate",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECDUPLICATE")
    RETURN
999 ERRORSEXITS("PETSC_VECDUPLICATE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECDUPLICATE
    
  !
  !================================================================================================================================
  !

!  !>Buffer routine to the PETSc VecGetArray routine.
!  SUBROUTINE PETSC_VECGETARRAY(X,ARRAY,err,error,*)

!    !Argument Variables
!    TYPE(PetscVecType), INTENT(INOUT), TARGET :: X !<The vector to get the array of
!    REAL(DP), POINTER :: ARRAY(:) !<On exit, a pointer to the array of the vector
!    INTEGER(INTG), INTENT(OUT) :: err !<The error code
!    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
!    !Local Variables

!    ENTERS("PETSC_VECGETARRAY",err,error,*999)

!    IF(ASSOCIATED(ARRAY)) THEN
!      CALL FlagError("Array is already associated",err,error,*999)
!    ELSE
!      CALL VecGetArray(X%vec,X%vec_DATA,X%vec_OFFSET,err)
!      IF(err/=0) THEN
!        IF(petscHandleError) THEN
!          CHKERRQ(err)
!        ENDIF
!        CALL FlagError("PETSc error in VecGetArray",err,error,*999)
!      ENDIF
!      ARRAY=>X%vec_DATA(X%vec_OFFSET:)
!    ENDIF
    
!    EXITS("PETSC_VECGETARRAY")
!    RETURN
!999 ERRORSEXITS("PETSC_VECGETARRAY",err,error)
!    EXITS("PETSC_VECGETARRAY")
!    RETURN 1
!  END SUBROUTINE PETSC_VECGETARRAY
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetArrayF90 routine.
  SUBROUTINE PETSC_VECGETARRAYF90(X,ARRAY,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT), TARGET :: X !<The vector to get the array of
    REAL(DP), POINTER :: ARRAY(:) !<On exit, a pointer to the array of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGETARRAYF90",err,error,*999)

    IF(ASSOCIATED(ARRAY)) THEN
      CALL FlagError("Array is already associated",err,error,*999)
    ELSE
      CALL VecGetArrayF90(X%vec,ARRAY,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in VecGetArrayF90",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("PETSC_VECGETARRAYF90")
    RETURN
999 ERRORSEXITS("PETSC_VECGETARRAYF90",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECGETARRAYF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetArrayReadF90 routine.
  SUBROUTINE Petsc_VecGetArrayReadF90(x,array,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT), TARGET :: x !<The vector to get the array of
    REAL(DP), POINTER :: array(:) !<On exit, a pointer to the array of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetArrayReadF90",err,error,*999)

    IF(ASSOCIATED(array)) THEN
      CALL FlagError("Array is already associated",err,error,*999)
    ELSE
      CALL VecGetArrayReadF90(x%vec,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in VecGetArrayReadF90",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_VecGetArrayReadF90")
    RETURN
999 ERRORSEXITS("Petsc_VecGetArrayReadF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetArrayReadF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetLocalSize routine.
  SUBROUTINE PETSC_VECGETLOCALSIZE(X,SIZE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to get the local size of
    INTEGER(INTG), INTENT(OUT) :: SIZE !<On exit, the local size of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGETLOCALSIZE",err,error,*999)

    CALL VecGetLocalSize(X%vec,SIZE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGetLocalSize",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGETLOCALSIZE")
    RETURN
999 ERRORSEXITS("PETSC_VECGETLOCALSIZE",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECGETLOCALSIZE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetOwnershipRange routine.
  SUBROUTINE PETSC_VECGETOWNERSHIPRANGE(X,LOW,HIGH,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to get the ownership range of 
    INTEGER(INTG), INTENT(OUT) :: LOW !<On exit, the low end of the range
    INTEGER(INTG), INTENT(OUT) :: HIGH !<On exit, the high end of the range
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGETOWNERSHIPRANGE",err,error,*999)

    CALL VecGetOwnershipRange(X%vec,LOW,HIGH,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
        ENDIF
      CALL FlagError("PETSc error in VecGetOwnershipRange",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGETOWNERSHIPRANGE")
    RETURN
999 ERRORSEXITS("PETSC_VECGETOWNERSHIPRANGE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGETOWNERSHIPRANGE

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDot routine.
  SUBROUTINE Petsc_VecDot(x,y,dotProduct,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(IN) :: x !<The vector x
    TYPE(PetscVecType), INTENT(IN) :: y !<The vector y
    REAL(DP), INTENT(OUT) :: dotProduct !<The dot product 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Petsc_VecDot",err,error,*999)
    
    CALL VecDot(x%vec,y%vec,dotProduct,err)

    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecDot",err,error,*999)
    ENDIF

    EXITS("Petsc_VecDot")
    RETURN
999 ERRORSEXITS("Petsc_VecDot",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecDot

    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecNorm routine.
  SUBROUTINE Petsc_VecNorm(x,normType,norm,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(IN) :: x !<The vector x to find the norm of
    NormType, INTENT(IN) :: normType !<The norm type
    REAL(DP), INTENT(OUT) :: norm !<On exit, the vector norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Petsc_VecNorm",err,error,*999)
    
    CALL VecNorm(x%vec,normType,norm,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecNorm",err,error,*999)
    ENDIF

    EXITS("Petsc_VecNorm")
    RETURN
999 ERRORSEXITS("Petsc_VecNorm",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecNorm

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetSize routine.
  SUBROUTINE PETSC_VECGETSIZE(X,SIZE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to get the size of
    INTEGER(INTG), INTENT(OUT) :: SIZE !<On exit, the size of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGETSIZE",err,error,*999)

    CALL VecGetSize(X%vec,SIZE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGetSize",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGETSIZE")
    RETURN
999 ERRORSEXITS("PETSC_VECGETSIZE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGETSIZE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetValues routine.
  SUBROUTINE PETSC_VECGETVALUES(X,N,INDICES,VALUES,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the values for
    INTEGER(INTG), INTENT(IN) :: N !<The number of indicies to get
    INTEGER(INTG), INTENT(IN) :: INDICES(*) !<The indices to get
    REAL(DP), INTENT(OUT) :: VALUES(*) !<On return, the values at the specified indicies
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGETVALUES",err,error,*999)

    CALL VecGetValues(X%vec,N,INDICES,VALUES,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGetValues",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGETVALUES")
    RETURN
999 ERRORSEXITS("PETSC_VECGETVALUES",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGETVALUES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostGetLocalForm routine.
  SUBROUTINE PETSC_VECGHOSTGETLOCALFORM(G,L,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: G !<The global form of the vector
    TYPE(PetscVecType), INTENT(INOUT) :: L !<On exit, the local form of the vector with ghosts
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGHOSTGETLOCALFORM",err,error,*999)

    CALL VecGhostGetLocalForm(G%vec,L%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostGetLocalForm",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGHOSTGETLOCALFORM")
    RETURN
999 ERRORSEXITS("PETSC_VECGHOSTGETLOCALFORM",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGHOSTGETLOCALFORM
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostRestoreLocalForm routine.
  SUBROUTINE PETSC_VECGHOSTRESTORELOCALFORM(G,L,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: G !<The global form of the vector
    TYPE(PetscVecType), INTENT(INOUT) :: L !<The local form of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGHOSTRESTORELOCALFORM",err,error,*999)

    CALL VecGhostRestoreLocalForm(G%vec,L%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostRestoreLocalForm",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGHOSTRESTORELOCALFORM")
    RETURN
999 ERRORSEXITS("PETSC_VECGHOSTRESTORELOCALFORM",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGHOSTRESTORELOCALFORM
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostUpdateBegin routine.
  SUBROUTINE PETSC_VECGHOSTUPDATEBEGIN(X,INSERT_MODE,SCATTER_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to begin the ghost update for
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    ScatterMode, INTENT(IN) :: SCATTER_MODE !<The scatter mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGHOSTUPDATEBEGIN",err,error,*999)

    CALL VecGhostUpdateBegin(X%vec,INSERT_MODE,SCATTER_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostUpdateBegin",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGHOSTUPDATEBEGIN")
    RETURN
999 ERRORSEXITS("PETSC_VECGHOSTUPDATEBEGIN",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGHOSTUPDATEBEGIN
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostUpdateEnd routine.
  SUBROUTINE PETSC_VECGHOSTUPDATEEND(X,INSERT_MODE,SCATTER_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to end the ghost update for
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    ScatterMode, INTENT(IN) :: SCATTER_MODE !<The scatter mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECGHOSTUPDATEEND",err,error,*999)

    CALL VecGhostUpdateEnd(X%vec,INSERT_MODE,SCATTER_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostUpdateEnd",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECGHOSTUPDATEEND")
    RETURN
999 ERRORSEXITS("PETSC_VECGHOSTUPDATEEND",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECGHOSTUPDATEEND
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecRestoreArrayF90 routine.
  SUBROUTINE PETSC_VECRESTOREARRAYF90(X,ARRAY,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to restore the array of
    REAL(DP), POINTER :: ARRAY(:) !<A pointer to the data to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECRESTOREARRAYF90",err,error,*999)

    CALL VecRestoreArrayF90(X%vec,ARRAY,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecRestoreArrayF90",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECRESTOREARRAYF90")
    RETURN
999 ERRORSEXITS("PETSC_VECRESTOREARRAYF90",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECRESTOREARRAYF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecRestoreArrayReadF90 routine.
  SUBROUTINE Petsc_VecRestoreArrayReadF90(x,array,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT), TARGET :: x !<The vector to restore the array for
    REAL(DP), POINTER :: array(:) !<A pointer to the array of the vector to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecRestoreArrayReadF90",err,error,*999)

    CALL VecRestoreArrayReadF90(x%vec,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecRestoreArrayReadF90",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecRestoreArrayReadF90")
    RETURN
999 ERRORSEXITS("Petsc_VecRestoreArrayReadF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecRestoreArrayReadF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecScale routine.
  SUBROUTINE PETSC_VECSCALE(X,ALPHA,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to scale
    REAL(DP), INTENT(IN) :: ALPHA !<The scaling factor
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSCALE",err,error,*999)

    CALL VecScale(X%vec,ALPHA,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecScale",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSCALE")
    RETURN
999 ERRORSEXITS("PETSC_VECSCALE",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECSCALE
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSet routine.
  SUBROUTINE PETSC_VECSET(X,VALUE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the value of
    REAL(DP), INTENT(IN) :: VALUE !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSET",err,error,*999)

    CALL VecSet(X%vec,VALUE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSet",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSET")
    RETURN
999 ERRORSEXITS("PETSC_VECSET",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECSET
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetFromOptions routine.
  SUBROUTINE PETSC_VECSETFROMOPTIONS(X,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSETFROMOPTIONS",err,error,*999)

    CALL VecSetFromOptions(X%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetFromOptions",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSETFROMOPTIONS")
    RETURN
999 ERRORSEXITS("PETSC_VECSETFROMOPTIONS",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECSETFROMOPTIONS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetLocalToGlobalMapping routine.
  SUBROUTINE PETSC_VECSETLOCALTOGLOBALMAPPING(X,CTX,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the local to global mapping for
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: CTX !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSETLOCALTOGLOBALMAPPING",err,error,*999)

    CALL VecSetLocalToGlobalMapping(X%vec,CTX%ISLOCALTOGLOBALMAPPING,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetLocalToGlobalMapping",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSETLOCALTOGLOBALMAPPING")
    RETURN
999 ERRORSEXITS("PETSC_VECSETLOCALTOGLOBALMAPPING",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECSETLOCALTOGLOBALMAPPING
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetValues routine.
  SUBROUTINE PETSC_VECSETVALUES(X,N,INDICES,VALUES,INSERT_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the values for
    INTEGER(INTG), INTENT(IN) :: N !<The number of indicies
    INTEGER(INTG), INTENT(IN) :: INDICES(*) !<The indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode, INTENT(IN) :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSETVALUES",err,error,*999)

    CALL VecSetValues(X%vec,N,INDICES,VALUES,INSERT_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetValues",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSETVALUES")
    RETURN
999 ERRORSEXITS("PETSC_VECSETVALUES",err,error)
    RETURN 1
  END SUBROUTINE PETSC_VECSETVALUES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetValuesLocal routine.
  SUBROUTINE PETSC_VECSETVALUESLOCAL(X,N,INDICES,VALUES,INSERT_MODE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the values of
    INTEGER(INTG), INTENT(IN) :: N !<The number of indices
    INTEGER(INTG), INTENT(IN) :: INDICES(*) !<The local indices
    REAL(DP), INTENT(IN) :: VALUES(*) !<The values to set
    InsertMode :: INSERT_MODE !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSETVALUESLOCAL",err,error,*999)

    CALL VecSetValuesLocal(X%vec,N,INDICES,VALUES,INSERT_MODE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetValuesLocal",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSETVALUESLOCAL")
    RETURN
999 ERRORSEXITS("PETSC_VECSETVALUESLOCAL",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECSETVALUESLOCAL
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetSizes routine.
  SUBROUTINE PETSC_VECSETSIZES(X,LOCAL_SIZE,GLOBAL_SIZE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to set the sizes of
    INTEGER(INTG), INTENT(IN) :: LOCAL_SIZE !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: GLOBAL_SIZE !<The number of global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECSETSIZES",err,error,*999)

    CALL VecSetSizes(X%vec,LOCAL_SIZE,GLOBAL_SIZE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetSizes",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECSETSIZES")
    RETURN
999 ERRORSEXITS("PETSC_VECSETSIZES",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECSETSIZES
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecView routine.
  SUBROUTINE PETSC_VECVIEW(X,V,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: X !<The vector to view
    PetscViewer, INTENT(IN) :: V !<The viewer
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("PETSC_VECVIEW",err,error,*999)

    CALL VecView(X%vec,V,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecView",err,error,*999)
    ENDIF
    
    EXITS("PETSC_VECVIEW")
    RETURN
999 ERRORSEXITS("PETSC_VECVIEW",err,error)
    RETURN 1
    
  END SUBROUTINE PETSC_VECVIEW


  !
  !================================================================================================================================
  !

END MODULE CMISS_PETSC
    
!>Buffer routine to the PETSc SNESSetJacobian routine for MatFDColoring contexts. The buffer is required because we want to
!>provide an interface so that we can pass a pointer to the solver for analytic Jacobian's. However, if we provided an interface
!>the Fortran's strong typing rules would not let us pass the matfdcoloring.
SUBROUTINE SNESSetJacobianBuffer(snes,A,B,jFunction,ctx,err)

  USE CmissPetscTypes
  USE Kinds

  IMPLICIT NONE
  
  !Argument Variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The snes to set the function for
  TYPE(PetscMatType), INTENT(INOUT) :: A !<The Jacobian matrix
  TYPE(PetscMatType), INTENT(INOUT) :: B !<The Jacobian preconditioning matrix
  EXTERNAL jFunction !<The external function to call
  TYPE(PetscMatFDColoringType) :: ctx !<The MatFDColoring data to pass to the function
  INTEGER(INTG), INTENT(OUT) :: err !<The error code

  CALL SNESSetJacobian(snes%snes,A%mat,B%mat,jFunction,ctx%matFDColoring,err)

END SUBROUTINE SNESSetJacobianBuffer
