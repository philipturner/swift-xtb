//
//  xTB_Calculator.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  /// Required.
  public var atomicNumbers: [UInt8]?
  
  /// Required.
  public var environment: xTB_Environment?
  
  /// Required.
  ///
  /// The default value is GFN2-xTB.
  public var hamiltonian: xTB_Hamiltonian = .tightBinding(2)
  
  /// Required. The net charge of the system.
  ///
  /// The default value is zero.
  public var netCharge: Float = .zero
  
  /// Required. The net spin of the system.
  ///
  /// The default value is zero.
  public var netSpin: Float = .zero
  
  /// Optional. The initial positions.
  ///
  /// When using GFN-FF, the positions are required to initialize force field
  /// parameters.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

/// Singlepoint calculator.
public class xTB_Calculator {
  // Keep a reference to the environment, so it is never deallocated while the
  // molecule is still in use.
  var environment: xTB_Environment
  
  var molecule: xTB_Molecule
  
  var calc: xtb_TCalculator?
  
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
}
