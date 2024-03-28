//
//  BLAS.swift
//
//
//  Created by Philip Turner on 3/28/24.
//

import Accelerate

// Make this library into a trampoline for Accelerate. Here are the symbols
// needed by xTB:
//
// This is a temporary solution for accelerating semiempirical simulations. The
// final solution should be full GPU acceleration of the remaining O(n^2)
// bottlenecks. That may require conversion of all FP64 operations to FP32,
// and full GPU offloading of the entire computation. It could be a while until
// that is practical. Therefore, only the linear algebra parts are accelerated
// for now, and it might not be from GPU acceleration.

/*

 00000000000c2bf0 T ___xtb_mctc_blas_level1_MOD_mctc_dasum
 00000000000c27b0 T ___xtb_mctc_blas_level1_MOD_mctc_daxpy
 00000000000c2390 T ___xtb_mctc_blas_level1_MOD_mctc_dcopy
 00000000000c1fd0 T ___xtb_mctc_blas_level1_MOD_mctc_ddot
 00000000000c1e10 T ___xtb_mctc_blas_level1_MOD_mctc_dnrm2
 00000000000c1934 T ___xtb_mctc_blas_level1_MOD_mctc_drot
 00000000000c16f0 T ___xtb_mctc_blas_level1_MOD_mctc_dscal
 00000000000c1250 T ___xtb_mctc_blas_level1_MOD_mctc_dswap
 00000000000c1030 T ___xtb_mctc_blas_level1_MOD_mctc_idamax
 00000000000c1140 T ___xtb_mctc_blas_level1_MOD_mctc_isamax
 00000000000c2cd0 T ___xtb_mctc_blas_level1_MOD_mctc_sasum
 00000000000c29d0 T ___xtb_mctc_blas_level1_MOD_mctc_saxpy
 00000000000c25a0 T ___xtb_mctc_blas_level1_MOD_mctc_scopy
 00000000000c21b0 T ___xtb_mctc_blas_level1_MOD_mctc_sdot
 00000000000c1ef0 T ___xtb_mctc_blas_level1_MOD_mctc_snrm2
 00000000000c1ba0 T ___xtb_mctc_blas_level1_MOD_mctc_srot
 00000000000c1814 T ___xtb_mctc_blas_level1_MOD_mctc_sscal
 00000000000c14a0 T ___xtb_mctc_blas_level1_MOD_mctc_sswap
 00000000000c5a50 T ___xtb_mctc_blas_level2_MOD_mctc_dgemv
 00000000000c52f0 T ___xtb_mctc_blas_level2_MOD_mctc_dger
 00000000000c4d50 T ___xtb_mctc_blas_level2_MOD_mctc_dspmv
 00000000000c4894 T ___xtb_mctc_blas_level2_MOD_mctc_dspr
 00000000000c4310 T ___xtb_mctc_blas_level2_MOD_mctc_dspr2
 00000000000c3c10 T ___xtb_mctc_blas_level2_MOD_mctc_dsymv
 00000000000c3530 T ___xtb_mctc_blas_level2_MOD_mctc_dsyr
 00000000000c2db0 T ___xtb_mctc_blas_level2_MOD_mctc_dsyr2
 00000000000c5dd4 T ___xtb_mctc_blas_level2_MOD_mctc_sgemv
 00000000000c56a0 T ___xtb_mctc_blas_level2_MOD_mctc_sger
 00000000000c5020 T ___xtb_mctc_blas_level2_MOD_mctc_sspmv
 00000000000c4af0 T ___xtb_mctc_blas_level2_MOD_mctc_sspr
 00000000000c45d4 T ___xtb_mctc_blas_level2_MOD_mctc_sspr2
 00000000000c3f90 T ___xtb_mctc_blas_level2_MOD_mctc_ssymv
 00000000000c38a0 T ___xtb_mctc_blas_level2_MOD_mctc_ssyr
 00000000000c3170 T ___xtb_mctc_blas_level2_MOD_mctc_ssyr2
 00000000000c94a4 T ___xtb_mctc_blas_level3_MOD_mctc_dgemm
 00000000000c89a4 T ___xtb_mctc_blas_level3_MOD_mctc_dsymm
 00000000000c74a4 T ___xtb_mctc_blas_level3_MOD_mctc_dsyr2k
 00000000000c8004 T ___xtb_mctc_blas_level3_MOD_mctc_dsyrk
 00000000000c6160 T ___xtb_mctc_blas_level3_MOD_mctc_dtrmm
 00000000000c6b04 T ___xtb_mctc_blas_level3_MOD_mctc_dtrsm
 00000000000c9a50 T ___xtb_mctc_blas_level3_MOD_mctc_sgemm
 00000000000c8f24 T ___xtb_mctc_blas_level3_MOD_mctc_ssymm
 00000000000c7a54 T ___xtb_mctc_blas_level3_MOD_mctc_ssyr2k
 00000000000c84d4 T ___xtb_mctc_blas_level3_MOD_mctc_ssyrk
 00000000000c6634 T ___xtb_mctc_blas_level3_MOD_mctc_strmm
 00000000000c6fd4 T ___xtb_mctc_blas_level3_MOD_mctc_strsm
 00000000000cb1f0 T ___xtb_mctc_blas_wrap1_MOD_mctc_dasum2
 00000000000caf30 T ___xtb_mctc_blas_wrap1_MOD_mctc_daxpy12
 00000000000cae80 T ___xtb_mctc_blas_wrap1_MOD_mctc_daxpy21
 00000000000cadd0 T ___xtb_mctc_blas_wrap1_MOD_mctc_daxpy22
 00000000000cab10 T ___xtb_mctc_blas_wrap1_MOD_mctc_dcopy12
 00000000000caa60 T ___xtb_mctc_blas_wrap1_MOD_mctc_dcopy21
 00000000000ca9b0 T ___xtb_mctc_blas_wrap1_MOD_mctc_dcopy22
 00000000000ca6f0 T ___xtb_mctc_blas_wrap1_MOD_mctc_ddot12
 00000000000ca640 T ___xtb_mctc_blas_wrap1_MOD_mctc_ddot21
 00000000000ca590 T ___xtb_mctc_blas_wrap1_MOD_mctc_ddot22
 00000000000ca4a0 T ___xtb_mctc_blas_wrap1_MOD_mctc_dnrm22
 00000000000ca340 T ___xtb_mctc_blas_wrap1_MOD_mctc_drot22
 00000000000ca250 T ___xtb_mctc_blas_wrap1_MOD_mctc_dscal2
 00000000000ca0f0 T ___xtb_mctc_blas_wrap1_MOD_mctc_dswap22
 00000000000ca000 T ___xtb_mctc_blas_wrap1_MOD_mctc_idamax2
 00000000000ca074 T ___xtb_mctc_blas_wrap1_MOD_mctc_isamax2
 00000000000cb264 T ___xtb_mctc_blas_wrap1_MOD_mctc_sasum2
 00000000000cb140 T ___xtb_mctc_blas_wrap1_MOD_mctc_saxpy12
 00000000000cb090 T ___xtb_mctc_blas_wrap1_MOD_mctc_saxpy21
 00000000000cafe0 T ___xtb_mctc_blas_wrap1_MOD_mctc_saxpy22
 00000000000cad20 T ___xtb_mctc_blas_wrap1_MOD_mctc_scopy12
 00000000000cac70 T ___xtb_mctc_blas_wrap1_MOD_mctc_scopy21
 00000000000cabc0 T ___xtb_mctc_blas_wrap1_MOD_mctc_scopy22
 00000000000ca900 T ___xtb_mctc_blas_wrap1_MOD_mctc_sdot12
 00000000000ca850 T ___xtb_mctc_blas_wrap1_MOD_mctc_sdot21
 00000000000ca7a0 T ___xtb_mctc_blas_wrap1_MOD_mctc_sdot22
 00000000000ca514 T ___xtb_mctc_blas_wrap1_MOD_mctc_snrm22
 00000000000ca3f0 T ___xtb_mctc_blas_wrap1_MOD_mctc_srot22
 00000000000ca2c4 T ___xtb_mctc_blas_wrap1_MOD_mctc_sscal2
 00000000000ca1a0 T ___xtb_mctc_blas_wrap1_MOD_mctc_sswap22
 00000000000cb474 T ___xtb_mctc_blas_wrap2_MOD_mctc_dgemv312
 00000000000cb2e0 T ___xtb_mctc_blas_wrap2_MOD_mctc_dgemv321
 00000000000cb7a4 T ___xtb_mctc_blas_wrap2_MOD_mctc_sgemv312
 00000000000cb610 T ___xtb_mctc_blas_wrap2_MOD_mctc_sgemv321
 00000000000cbbd4 T ___xtb_mctc_blas_wrap3_MOD_mctc_dgemm233
 00000000000cbdd0 T ___xtb_mctc_blas_wrap3_MOD_mctc_dgemm323
 00000000000cb940 T ___xtb_mctc_blas_wrap3_MOD_mctc_dgemm332
 00000000000cc264 T ___xtb_mctc_blas_wrap3_MOD_mctc_sgemm233
 00000000000cc3e0 T ___xtb_mctc_blas_wrap3_MOD_mctc_sgemm323
 00000000000cbfd0 T ___xtb_mctc_blas_wrap3_MOD_mctc_sgemm332

 00000000000cc5e0 T ___xtb_mctc_lapack_eigensolve_MOD___copy_xtb_mctc_lapack_eigensolve_Teigensolver
 000000000034c2a8 S ___xtb_mctc_lapack_eigensolve_MOD___def_init_xtb_mctc_lapack_eigensolve_Teigensolver
 00000000000cc770 T ___xtb_mctc_lapack_eigensolve_MOD___final_xtb_mctc_lapack_eigensolve_Teigensolver
 000000000049e518 D ___xtb_mctc_lapack_eigensolve_MOD___vtab_xtb_mctc_lapack_eigensolve_Teigensolver
 00000000000cd3d0 T ___xtb_mctc_lapack_eigensolve_MOD_initdeigensolver
 00000000000cd764 T ___xtb_mctc_lapack_eigensolve_MOD_initseigensolver
 00000000000cc950 T ___xtb_mctc_lapack_eigensolve_MOD_mctc_dsygvd
 00000000000cce90 T ___xtb_mctc_lapack_eigensolve_MOD_mctc_ssygvd
 00000000000cdb00 T ___xtb_mctc_lapack_gst_MOD_mctc_dspgst
 00000000000ce0c0 T ___xtb_mctc_lapack_gst_MOD_mctc_dsygst
 00000000000cdde0 T ___xtb_mctc_lapack_gst_MOD_mctc_sspgst
 00000000000ce540 T ___xtb_mctc_lapack_gst_MOD_mctc_ssygst
 00000000000d0620 T ___xtb_mctc_lapack_trf_MOD_mctc_dgetrf
 00000000000ced20 T ___xtb_mctc_lapack_trf_MOD_mctc_dpotrf
 00000000000ce9c0 T ___xtb_mctc_lapack_trf_MOD_mctc_dpptrf
 00000000000cf220 T ___xtb_mctc_lapack_trf_MOD_mctc_dsptrf
 00000000000cf800 T ___xtb_mctc_lapack_trf_MOD_mctc_dsytrf
 00000000000d09b0 T ___xtb_mctc_lapack_trf_MOD_mctc_sgetrf
 00000000000cefa0 T ___xtb_mctc_lapack_trf_MOD_mctc_spotrf
 00000000000ceb70 T ___xtb_mctc_lapack_trf_MOD_mctc_spptrf
 00000000000cf510 T ___xtb_mctc_lapack_trf_MOD_mctc_ssptrf
 00000000000cff10 T ___xtb_mctc_lapack_trf_MOD_mctc_ssytrf
 00000000000d2324 T ___xtb_mctc_lapack_tri_MOD_mctc_dgetri
 00000000000d10a0 T ___xtb_mctc_lapack_tri_MOD_mctc_dpotri
 00000000000d0d40 T ___xtb_mctc_lapack_tri_MOD_mctc_dpptri
 00000000000d15a0 T ___xtb_mctc_lapack_tri_MOD_mctc_dsptri
 00000000000d1b84 T ___xtb_mctc_lapack_tri_MOD_mctc_dsytri
 00000000000d2980 T ___xtb_mctc_lapack_tri_MOD_mctc_sgetri
 00000000000d1320 T ___xtb_mctc_lapack_tri_MOD_mctc_spotri
 00000000000d0ef0 T ___xtb_mctc_lapack_tri_MOD_mctc_spptri
 00000000000d1894 T ___xtb_mctc_lapack_tri_MOD_mctc_ssptri
 00000000000d1f54 T ___xtb_mctc_lapack_tri_MOD_mctc_ssytri
 00000000000d5220 T ___xtb_mctc_lapack_trs_MOD_mctc_dgetrs
 00000000000d3770 T ___xtb_mctc_lapack_trs_MOD_mctc_dpotrs
 00000000000d2fe0 T ___xtb_mctc_lapack_trs_MOD_mctc_dpptrs
 00000000000d4070 T ___xtb_mctc_lapack_trs_MOD_mctc_dsptrs
 00000000000d48d4 T ___xtb_mctc_lapack_trs_MOD_mctc_dsytrs
 00000000000d56c0 T ___xtb_mctc_lapack_trs_MOD_mctc_sgetrs
 00000000000d3bf0 T ___xtb_mctc_lapack_trs_MOD_mctc_spotrs
 00000000000d33b0 T ___xtb_mctc_lapack_trs_MOD_mctc_spptrs
 00000000000d44b4 T ___xtb_mctc_lapack_trs_MOD_mctc_ssptrs
 00000000000d4d80 T ___xtb_mctc_lapack_trs_MOD_mctc_ssytrs
 00000000000d6d90 T ___xtb_mctc_lapack_wrap_MOD_mctc_dgetrs1
 00000000000d6c30 T ___xtb_mctc_lapack_wrap_MOD_mctc_dgetrs3
 00000000000d5fc4 T ___xtb_mctc_lapack_wrap_MOD_mctc_dpotrs1
 00000000000d5ec0 T ___xtb_mctc_lapack_wrap_MOD_mctc_dpotrs3
 00000000000d5c50 T ___xtb_mctc_lapack_wrap_MOD_mctc_dpptrs1
 00000000000d5b60 T ___xtb_mctc_lapack_wrap_MOD_mctc_dpptrs3
 00000000000d63c0 T ___xtb_mctc_lapack_wrap_MOD_mctc_dsptrs1
 00000000000d6280 T ___xtb_mctc_lapack_wrap_MOD_mctc_dsptrs3
 00000000000d6880 T ___xtb_mctc_lapack_wrap_MOD_mctc_dsytrs1
 00000000000d6720 T ___xtb_mctc_lapack_wrap_MOD_mctc_dsytrs3
 00000000000d7020 T ___xtb_mctc_lapack_wrap_MOD_mctc_sgetrs1
 00000000000d6ec0 T ___xtb_mctc_lapack_wrap_MOD_mctc_sgetrs3
 00000000000d61a4 T ___xtb_mctc_lapack_wrap_MOD_mctc_spotrs1
 00000000000d60a0 T ___xtb_mctc_lapack_wrap_MOD_mctc_spotrs3
 00000000000d5e00 T ___xtb_mctc_lapack_wrap_MOD_mctc_spptrs1
 00000000000d5d10 T ___xtb_mctc_lapack_wrap_MOD_mctc_spptrs3
 00000000000d6610 T ___xtb_mctc_lapack_wrap_MOD_mctc_ssptrs1
 00000000000d64d0 T ___xtb_mctc_lapack_wrap_MOD_mctc_ssptrs3
 00000000000d6b10 T ___xtb_mctc_lapack_wrap_MOD_mctc_ssytrs1
 00000000000d69b0 T ___xtb_mctc_lapack_wrap_MOD_mctc_ssytrs3
 
                  U _dgemm_
                  U _dgemv_
                  U _dger_
                  U _dgetrf_
                  U _dgetri_
                  U _dgetrs_
                  U _dnrm2_
                  U _dpotrf_
                  U _dpotri_
                  U _dpotrs_
                  U _dpptrf_
                  U _dpptri_
                  U _dpptrs_
                  U _drot_
                  U _dscal_
                  U _dspev_
                  U _dspgst_
                  U _dspmv_
                  U _dspr2_
                  U _dspr_
                  U _dsptrf_
                  U _dsptri_
                  U _dsptrs_
                  U _dswap_
                  U _dsyev_
                  U _dsyevd_
                  U _dsygst_
                  U _dsygvd_
                  U _dsymm_
                  U _dsymv_
                  U _dsyr2_
                  U _dsyr2k_
                  U _dsyr_
                  U _dsyrk_
                  U _dsysv_
                  U _dsytrf_
                  U _dsytri_
                  U _dsytrs_
                  U _dtrmm_
                  U _dtrsm_

                  U _sasum_
                  U _saxpy_
                  U _scopy_
                  U _sdot_
                  U _sgemm_
                  U _sgemv_
                  U _sger_
                  U _sgetrf_
                  U _sgetri_
                  U _sgetrs_
                  U _snrm2_
                  U _spotrf_
                  U _spotri_
                  U _spotrs_
                  U _spptrf_
                  U _spptri_
                  U _spptrs_
                  U _sqrt
                  U _srot_
                  U _sscal_
                  U _sscanf
                  U _sspevx_
                  U _sspgst_
                  U _sspmv_
                  U _sspr2_
                  U _sspr_
                  U _ssptrf_
                  U _ssptri_
                  U _ssptrs_
                  U _sswap_
                  U _ssyev_
                  U _ssyevd_
                  U _ssyevx_
                  U _ssygst_
                  U _ssygvd_
                  U _ssymm_
                  U _ssymv_
                  U _ssyr2_
                  U _ssyr2k_
                  U _ssyr_
                  U _ssyrk_
                  U _ssytrf_
                  U _ssytri_
                  U _ssytrs_
                  U _strmm_
                  U _strsm_
 */

// TODO: Create a Swift script that can translate editor placeholders into
// overridden BLAS functions. Store the script somewhere in the Mechanosynthesis
// docs, perhaps in a new folder about accelerating xTB.

//@_cdecl("sgemm_")
//func sgemm_(_ TRANSA: UnsafePointer<CChar>, _ TRANSB: UnsafePointer<CChar>, _ M: UnsafePointer<__LAPACK_int>, _ N: UnsafePointer<__LAPACK_int>, _ K: UnsafePointer<__LAPACK_int>, _ ALPHA: UnsafePointer<Float>, _ A: UnsafePointer<Float>?, _ LDA: UnsafePointer<__LAPACK_int>, _ B: UnsafePointer<Float>?, _ LDB: UnsafePointer<__LAPACK_int>, _ BETA: UnsafePointer<Float>, _ C: UnsafeMutablePointer<Float>?, _ LDC: UnsafePointer<__LAPACK_int>) {
//  Accelerate.sgemm_(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
//}
