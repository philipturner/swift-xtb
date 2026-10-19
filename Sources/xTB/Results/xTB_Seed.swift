//
//  xTB_Seed.swift
//  swift-xtb
//
//  Created by Philip Turner on 10/18/26.
//

/// A converged SCF solution, used to warm start a future singlepoint.
///
/// Create a seed with ``xTB_Calculator/seed``, then apply it with
/// ``xTB_Calculator/withSeed(_:_:)``. A seed can be applied any number of
/// times, but only to the calculator it came from.
public class xTB_Seed {
  public enum Error: Swift.Error, Equatable, CustomStringConvertible {
    case differentCalculator
    case unused
    
    public var description: String {
      switch self {
      case .differentCalculator:
        return "The seed was created by a different calculator."
      case .unused:
        return "No singlepoint was computed, so the seed was never used."
      }
    }
  }
  
  weak var calculator: xTB_Calculator?
  
  let results: xTB_Results
  
  init(calculator: xTB_Calculator, results: xTB_Results) {
    self.calculator = calculator
    self.results = xTB_Results(copying: results)
  }
}

extension xTB_Calculator {
  /// Captures the current SCF solution.
  ///
  /// Returns `nil` if there is no up-to-date singlepoint. That happens before
  /// the first property (energy, forces, etc.) is evaluated, and after any
  /// change to the calculator, until a property is evaluated again.
  public var seed: xTB_Seed? {
    guard let results else {
      return nil
    }
    return xTB_Seed(calculator: self, results: results)
  }
  
  /// Warm starts every singlepoint computed inside `body`.
  ///
  /// Each singlepoint starts from the seed, not from the previous singlepoint
  /// inside `body`. After `body` returns, singlepoints fall back to a cold
  /// start.
  ///
  /// - Throws: ``xTB_Seed/Error/differentCalculator`` before running `body`,
  ///   if the seed came from another calculator.
  ///   ``xTB_Seed/Error/unused`` after running `body`, if no singlepoint was
  ///   computed. For example, the body read a property that was already
  ///   cached.
  public func withSeed<R>(
    _ seed: xTB_Seed,
    _ body: () throws -> R
  ) throws -> R {
    guard seed.calculator === self else {
      throw xTB_Seed.Error.differentCalculator
    }
    guard pendingSeed == nil else {
      fatalError("Calls to withSeed cannot be nested.")
    }
    
    pendingSeed = seed
    pendingSeedUsed = false
    defer { pendingSeed = nil }
    let output = try body()
    
    guard pendingSeedUsed else {
      throw xTB_Seed.Error.unused
    }
    return output
  }
}
