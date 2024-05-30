//
//  xTB_Calculator.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  public var atomicNumbers: [UInt8]?
  
  public var environment: xTB_Environment?
  
  public var netCharge: Float = .zero
  
  public var netSpin: Float = .zero
  
  public init() {
    
  }
}

/// Singlepoint calculator.
public class xTB_Calculator {
  var calc: xtb_TCalculator?
  
  var molecule: xTB_Molecule
  
  var environment: xTB_Environment { molecule.environment }
  
  /// Create new calculator object.
  public init(descriptor: xTB_CalculatorDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let environment = descriptor.environment else {
      fatalError("Descriptor was incomplete.")
    }
    
    // Create the backing molecule.
    var moleculeDesc = xTB_MoleculeDescriptor()
    moleculeDesc.atomicNumbers = atomicNumbers
    moleculeDesc.environment = environment
    moleculeDesc.netCharge = descriptor.netCharge
    moleculeDesc.netSpin = descriptor.netSpin
    molecule = xTB_Molecule(descriptor: moleculeDesc)
    
    // Create the backing calculator.
    calc = xtb_newCalculator()
    guard calc != nil else {
      fatalError("Could not create new xTB_Calculator.")
    }
    
    // Initialize GFN2-xTB.
    xtb_loadGFN2xTB(
      environment.env, molecule.mol, calc, nil)
  }
  
  /// Delete calculator object.
  deinit {
    xtb_delCalculator(&calc)
  }
  
  /// Update coordinates (in nm).
  public func setPositions(_ positions: [SIMD3<Float>]) {
    molecule.setPositions(positions)
  }
  
  // TODO: Find an ergonomic API for point charges. Perhaps as a separate
  // class, whose gradient can be queried because it stores a strong reference
  // to the calculator object.
}
