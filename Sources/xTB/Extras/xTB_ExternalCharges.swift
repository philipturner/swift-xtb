//
//  xTB_ExternalCharges.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

/// External charge potential.
public struct xTB_ExternalCharge {
  /// Partial charge in units of proton charge.
  public var charge: Float = .zero
  
  /// Chemical hardness (in atomic units).
  public var chemicalHardness: UInt32 = .zero
  
  /// Coordinates (in nm).
  public var position: SIMD3<Float> = .zero
  
  public init() {
    
  }
}

extension xTB_Calculator {
  public func setExternalCharges(_ externalCharges: [xTB_ExternalCharge]?) {
    
  }
  
  func releaseExternalCharges() {
    
  }
}

/// Add a external charge potential to calculator (only supported in GFN1/2-xTB)
let xtb_setExternalCharges: @convention(c) (
  xtb_TEnvironment?,
  xtb_TCalculator?,
  UnsafeMutablePointer<Int32>?, // n
  UnsafeMutablePointer<Int32>?, // numbers [n]
  UnsafeMutablePointer<Double>?, // charges [n]
  UnsafeMutablePointer<Double> // positions [n][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setExternalCharges")
