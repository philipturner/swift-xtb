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
extension xTB_Environment {
  nonisolated(unsafe)
  private static var _threadCount: Int?

  /// The number of threads for xTB's parallel regions (integrals, gradients,
  /// ...), or nil to leave the OpenMP runtime's own setting: OMP_NUM_THREADS
  /// as it was when the program launched, else `maximumThreadCount`.
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
