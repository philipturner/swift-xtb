//
//  xTB_Calculator.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/29/24.
//

import C_xTB

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. The parameterized method for evaluating forces.
  ///
  /// The default value is GFN2-xTB.
  public var hamiltonian: xTB_Hamiltonian = .tightBinding
  
  /// Required. The net charge of the system.
  ///
  /// The default value is zero.
  public var netCharge: Float = .zero
  
  /// Required. The net spin of the system.
  ///
  /// The default value is zero.
  public var netSpin: Float = .zero
  
  /// Optional. The position of each atom's nucleus (in nanometers).
  ///
  /// When using GFN-FF, the positions are needed to initialize force field
  /// parameters. When using GFN2-xTB, positions can be specified after
  /// initialization.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

/// Singlepoint calculator.
public class xTB_Calculator {
  /// The parameterized method for evaluating forces.
  public let hamiltonian: xTB_Hamiltonian
  
  var tCalculator: xtb_TCalculator!
  var tMolecule: xtb_TMolecule!
  
  var storage: xTB_CalculatorStorage
  var results: xTB_Results!
  
  public init(descriptor: xTB_CalculatorDescriptor) {
    self.hamiltonian = descriptor.hamiltonian
    
    let molecule = xTB_Molecule(descriptor: descriptor)
    let orbitals = xTB_Orbitals(descriptor: descriptor)
    
    // Create the 'TCalculator'.
    self.tCalculator = xTB_Calculator.createObject()
    
    // Create the 'TMolecule'.
    self.tMolecule = xTB_Molecule.createObject(molecule)
    
    // Assign ownership of 'molecule' and 'orbitals' to 'storage'.
    storage = xTB_CalculatorStorage()
    storage.molecule = molecule
    storage.orbitals = orbitals
    storage.molecule.calculator = self
    storage.orbitals.calculator = self
    
    // Load the parameters.
    switch descriptor.hamiltonian {
    case .forceField:
      xtb_loadGFNFF(
        xTB_Environment.tEnvironment, tMolecule, tCalculator, nil)
    case .tightBinding:
      xtb_loadGFN2xTB(
        xTB_Environment.tEnvironment, tMolecule, tCalculator, nil)
    }
  }
  
  deinit {
    xtb_delMolecule(&tMolecule)
    xtb_delCalculator(&tCalculator)
  }
  
  /// Create the reference-counted object from the C API.
  static func createObject() -> xtb_TCalculator {
    let calc = xtb_newCalculator()
    guard let calc else {
      fatalError("Could not create new xTB_Calculator.")
    }
    return calc
  }
}
