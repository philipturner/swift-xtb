//
//  xTB_ExternalCharges.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

/// External charge potential.
public struct xTB_ExternalCharges {
  unowned var calculator: xTB_Calculator?
  
  /// The atomic number to match the chemical hardness to.
  ///
  /// A good default value is 7 (nitrogen).
  public var atomicNumbers: [UInt8] = []
  
  /// Partial charge in units of proton charge.
  public var charges: [Float] = []
  
  /// The position of each point charge (in nm).
  public var positions: [SIMD3<Float>] = []
}

extension xTB_ExternalCharges {
  func update() {
    print("Updating external charges.")
    
    // Erase the previous external potential.
    xtb_releaseExternalCharges(
      xTB_Environment._environment,
      calculator!._calculator)
    
    // Check that each array has the same length.
    guard atomicNumbers.count == charges.count,
          atomicNumbers.count == positions.count else {
      fatalError("Arrays did not have the same length.")
    }
    if atomicNumbers.count == 0 {
      // Bypass an error when there are zero charges.
      return
    }
    
    // Set the external charge potential.
    var n = Int32(atomicNumbers.count)
    var numbers = atomicNumbers.map(Int32.init)
    var charges = self.charges.map(Double.init)
    var positions64 = xTB_Molecule.convertPositions(positions)
    
    // Initialize the external potential.
    xtb_setExternalCharges(
      xTB_Environment._environment,
      calculator!._calculator,
      &n,
      &numbers,
      &charges,
      &positions64)
  }
}
