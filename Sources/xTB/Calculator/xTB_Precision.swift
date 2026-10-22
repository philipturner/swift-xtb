//
//  xTB_Precision.swift
//  swift-xtb
//
//  Created by Philip Turner on 10/21/26.
//

/// The precision of the eigensolves in the self-consistent charge (SCC)
/// iterations of GFN2-xTB.
public enum xTB_Precision: Equatable, Sendable {
  /// Every eigensolve in FP64 (double precision), as upstream xTB does.
  case fp64

  /// FP32 (single precision) eigensolves while the root-mean-square change of
  /// the charges exceeds `switchThreshold`, then FP64 until convergence.
  ///
  /// The final iterations always run in FP64, so a converged result is an
  /// FP64 solution. It differs from the `.fp64` result by about as much as
  /// the SCC's own convergence noise (production accuracy: within 2e-7
  /// Hartree in energy over the test set). If the FP32 attempt hasn't
  /// converged after `maximumIterations` iterations, xTB reruns the SCC from
  /// the same start in FP64, which reproduces the `.fp64` result exactly.
  ///
  /// `refined`: each FP32 eigensolve is corrected in FP64, which keeps the
  /// results within about 1e-12 Hartree of `.fp64`. It is faster than `.fp64`
  /// for one process, but not for one process per core.
  ///
  /// Measured on an M1 Max (n = 412-776 orbitals): 12-20% faster per
  /// singlepoint than `.fp64`, 15-25% more throughput with 8 processes x 1
  /// thread.
  ///
  /// > Note: Requires an xTB library with the SCC precision API
  /// (``xTB_Library/supportsFP32``). Not available for GFN-FF.
  case fp32(
    switchThreshold: Double = 1e-3,
    maximumIterations: Int = 160,
    refined: Bool = false)
}

/// What the SCC of a singlepoint actually ran, as reported by xTB.
public struct xTB_PrecisionReport: Equatable, Sendable {
  /// The precision xTB used. Always equal to the calculator's
  /// ``xTB_Calculator/precision``: a mismatch is a fatal error.
  public var precision: xTB_Precision

  /// Iterations with an FP32 eigensolve (including those of an FP32 attempt
  /// that was then rerun in FP64).
  public var fp32Iterations: Int

  /// Whether the FP32 attempt failed to converge and xTB reran the SCC from
  /// the same start in FP64.
  public var restartedInFP64: Bool

  /// SCC iterations of the run that produced the results.
  public var iterations: Int
}

extension xTB_Calculator {
  /// Precision of the eigensolves in the SCC iterations.
  ///
  /// The default value is `.fp64`.
  ///
  /// > Note: `.fp32` requires an xTB library with the SCC precision API
  /// (``xTB_Library/supportsFP32``); requesting it from another library is a
  /// fatal error. Not available for GFN-FF.
  public var precision: xTB_Precision {
    get {
      storage.precision
    }
    set {
      applyPrecision(newValue)
      storage.precision = newValue
      invalidateSinglepoint()
    }
  }

  /// What the SCC of the current singlepoint ran: precision, FP32 iterations,
  /// FP64 restart.
  ///
  /// `nil` for GFN-FF, for xTB libraries without the SCC precision API, and
  /// when the singlepoint failed before its SCC.
  public var precisionReport: xTB_PrecisionReport? {
    requestSinglepoint()
    return results.precisionReport
  }

  /// Validate the precision and set it on the C API calculator.
  func applyPrecision(_ precision: xTB_Precision) {
    var threshold: Double = .zero
    var maximumIterations: Int = 160
    var refined: Bool = false
    switch precision {
    case .fp64:
      // Libraries without the API, and GFN-FF, always run FP64.
      guard xtb_supportsSCCPrecision, hamiltonian == .tightBinding else {
        return
      }
    case .fp32(let switchThreshold, let fp32MaximumIterations, let fp32Refined):
      guard xtb_supportsSCCPrecision else {
        fatalError("""
          FP32 eigensolves were requested, but the xTB library has no SCC \
          precision API (xtb_setSCCPrecision).
          """)
      }
      guard hamiltonian == .tightBinding else {
        fatalError("FP32 eigensolves are not available for GFN-FF.")
      }
      guard switchThreshold > 0, switchThreshold <= 1e3 else {
        fatalError("""
          FP32 switch threshold must be in (0, 1000], got \(switchThreshold).
          """)
      }
      guard fp32MaximumIterations >= 1,
            fp32MaximumIterations <= Int(Int32.max) else {
        fatalError("""
          FP32 maximum iterations must be positive, got \
          \(fp32MaximumIterations).
          """)
      }
      threshold = switchThreshold
      maximumIterations = fp32MaximumIterations
      refined = fp32Refined
    }

    // The values were validated above, so the C API should accept them. An
    // error that wasn't there before the call is fatal.
    let statusBefore = xTB_Environment.status
    xtb_setSCCPrecision!(
      xTB_Environment.tEnvironment,
      tCalculator,
      threshold,
      Int32(maximumIterations),
      refined ? 1 : 0)
    if statusBefore == 0, xTB_Environment.status != 0 {
      fatalError("""
        xtb_setSCCPrecision rejected \(precision): \
        \(xTB_Environment.flushErrorStack())
        """)
    }
  }
}

extension xTB_Results {
  /// Read what the SCC ran and check it against the calculator's precision.
  func checkPrecision() {
    guard xtb_supportsSCCPrecision,
          calculator.hamiltonian == .tightBinding else {
      return
    }

    var threshold: Double = .zero
    var maximumIterations: Int32 = .zero
    var refined: Int32 = .zero
    var fp32Iterations: Int32 = .zero
    var restarted: Int32 = .zero
    var iterations: Int32 = .zero
    xtb_getSCCPrecision!(
      xTB_Environment.tEnvironment, tResults,
      &threshold, &maximumIterations, &refined,
      &fp32Iterations, &restarted, &iterations)

    // No SCC ran: the singlepoint failed before it, and xTB reported an error.
    guard iterations > 0 else {
      return
    }

    let reported: xTB_Precision
    if threshold == 0 {
      reported = .fp64
    } else {
      reported = .fp32(
        switchThreshold: threshold,
        maximumIterations: Int(maximumIterations),
        refined: refined == 1)
    }
    guard reported == calculator.precision else {
      fatalError("""
        The calculator's precision is \(calculator.precision), but xTB ran \
        the SCC with \(reported).
        """)
    }
    precisionReport = xTB_PrecisionReport(
      precision: reported,
      fp32Iterations: Int(fp32Iterations),
      restartedInFP64: restarted == 1,
      iterations: Int(iterations))
  }
}
