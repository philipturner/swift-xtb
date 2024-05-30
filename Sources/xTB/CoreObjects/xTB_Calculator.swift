//
//  xTB_Calculator.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. The backing calculation environment.
  public var environment: xTB_Environment?
  
  /// Required. The parametrized method for evaluating forces.
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
  var pointer: xtb_TCalculator
  
  // The backing calculation environment.
  var environment: xTB_Environment
  
  // The backing molecule.
  var molecule: xTB_Molecule
  
  /// Create new calculator object.
  public init(descriptor: xTB_CalculatorDescriptor) {
    // Retain the backing calculation environment.
    guard let environment = descriptor.environment else {
      fatalError("Descriptor was incomplete.")
    }
    self.environment = environment
    
    // Create the backing calculator.
    guard let calc = xtb_newCalculator() else {
      fatalError("Could not create new xTB_Calculator.")
    }
    self.pointer = calc
    
    // Create the backing molecule.
    self.molecule = xTB_Molecule(descriptor: descriptor)
    
    // Initialize the positions.
    if let positions = descriptor.positions {
      setPositions(positions)
    }
    
    // Initialize the parameters.
    loadHamiltonian(descriptor.hamiltonian)
  }
  
  /// Delete calculator object.
  deinit {
    xtb_delCalculator(&pointer)
  }
}
