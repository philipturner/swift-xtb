//
//  xTB_Environment+OpenMP.swift
//  swift-xtb
//
//  Created by Philip Turner on 10/21/26.
//

import Foundation

// xTB parallelizes with OpenMP (GCC's libgomp in the default build). The
// runtime reads OMP_NUM_THREADS once, when it loads with the xTB library;
// changing the environment variable later has no effect. The count can then be
// changed only with omp_set_num_threads, which applies to the calling thread
// alone: a new Thread or a Dispatch worker starts from the load-time value.

/*
 An unresolved failure on Windows (October 2026)
 
 While this file was being written, OpenMPTests.testThreadCount failed on a
 Windows PC (x86_64-unknown-windows-msvc, MSYS2 GCC's libgomp-1.dll, 4-core
 Intel from 2014) for several rounds. Every team stayed at the runtime's
 default (4) after setting `threadCount` and running a calculation, on the
 calling thread, a new Thread and a Dispatch worker alike. Step-by-step checks
 (now testRuntimeLookup) narrowed it to applyThreadCount itself:
 
 - Exactly one OpenMP runtime was loaded: libgomp-1.dll in the package folder.
 - omp_set_num_threads called directly from the test worked (team 2).
 - xTB_Library.isLoaded was true, and loadOpenMPSymbol returned the same
   omp_set_num_threads address the test called successfully.
 - Yet applyThreadCount alone, on a fresh thread, left the team at 4.
 - xTB doesn't reset the count: the library imports omp_get_max_threads only,
   never omp_set_num_threads, so the singlepoint couldn't undo it.
 
 Deleting .build and rebuilding, many times, didn't help. What fixed it was
 purging the package folder and cloning it fresh: at the same commit as a
 failing run, every test passed, with no code change. A fresh clone also
 resets what deleting .build leaves alone (the DLLs in the package folder:
 xtb.dll, libgomp-1.dll and the other GCC runtime DLLs; Package.resolved;
 .swiftpm; local edits), so the cause was probably there, but which one was
 never found: the old folder was gone before it could be compared. Around the
 same time, one run also aborted with "libgomp: Thread creation failed:
 Resource temporarily unavailable" (test.sh sets OMP_STACKSIZE=2G); that
 didn't recur either.
 
 If the tests fail like this again on Windows: first check where the DLLs in
 the package folder came from and whether a fresh clone also fails, then read
 testRuntimeLookup's assertion messages, which name the first step that broke.
 */
extension xTB_Environment {
  nonisolated(unsafe)
  private static var _threadCount: Int?

  /// The number of threads for xTB's parallel regions (integrals, gradients,
  /// ...), or nil to leave the OpenMP runtime's own setting: OMP_NUM_THREADS
  /// as it was when the program launched, else `maximumThreadCount`.
  ///
  /// With the xTB library this branch builds (philipturner/xtb's
  /// 'diagonalization'), the count also covers the eigensolve in every SCC
  /// iteration, from 256 orbitals up: the divide and conquer uses all the
  /// threads, the reduction to tridiagonal form at most one cluster of the
  /// fastest cores (4 on an M1 Max), the back-transformation two. On an M1
  /// Max, a singlepoint was no faster with 8 threads than with 4.
  ///
  /// Accelerate's own threads (BLAS and LAPACK calls) follow the environment
  /// variable VECLIB_MAXIMUM_THREADS, which must be set before the program
  /// starts; swift-xtb doesn't change it. It barely matters: in every case
  /// tested (M1 Max, 8 OpenMP threads, n = 412-1204, FP64 and FP32), a
  /// singlepoint took within about 3% of the same time with 1 or 8 Accelerate
  /// threads, and Accelerate's divide and conquer alone changed by at most 4%.
  ///
  /// When set, swift-xtb sets the count on the calling thread before every
  /// singlepoint and parameter load, so calculations started on any thread
  /// use it. Setting it never loads the xTB library. Setting nil again stops
  /// the overrides, but threads that already ran a calculation keep the count
  /// they had last.
  public static var threadCount: Int? {
    get { _threadCount }
    set {
      if let newValue {
        precondition(newValue > 0, "The thread count must be positive.")
      }
      _threadCount = newValue
    }
  }

  /// The nominal number of cores (on Apple silicon, performance and efficiency
  /// cores together): the OpenMP runtime's thread count when OMP_NUM_THREADS
  /// isn't set. Not necessarily the fastest choice; set `threadCount` for
  /// that. Never loads the xTB library.
  public static var maximumThreadCount: Int {
    ProcessInfo.processInfo.activeProcessorCount
  }

  /// Sets the calling thread's OpenMP thread count to `threadCount`, if it is
  /// set. Call only with the xTB library loaded.
  static func applyThreadCount() {
    typealias Function = @convention(c) (Int32) -> Void
    guard let threadCount = _threadCount,
          let symbol = xTB_Library.loadOpenMPSymbol("omp_set_num_threads") else {
      return
    }
    unsafeBitCast(symbol, to: Function.self)(Int32(threadCount))
  }
}
