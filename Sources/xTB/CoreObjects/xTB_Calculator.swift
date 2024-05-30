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
  
  /// Optional. The position of each atom's nucleus (in nm).
  ///
  /// When using GFN-FF, the positions are required to initialize force field
  /// parameters. When using tight binding, positions can be specified
  /// after initialization.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

/// Singlepoint calculator.
public class xTB_Calculator {
  var pointer: xtb_TCalculator
  var environment: xTB_Environment
  public let molecule: xTB_Molecule
  public let orbitals: xTB_Orbitals
  var results: xTB_Results?
  
  var state = State()
  var updateRecord = UpdateRecord()
  
  /// Create new calculator object.
  public init(descriptor: xTB_CalculatorDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let environment = descriptor.environment else {
      fatalError("Descriptor was incomplete.")
    }
    self.environment = environment
    
    // Create the calculator.
    guard let calc = xtb_newCalculator() else {
      fatalError("Could not create new xTB_Calculator.")
    }
    self.pointer = calc
    
    // Create the initial positions.
    if let positions = descriptor.positions {
      state.positions = positions
    } else {
      state.positions = xTB_Molecule
        .createInitialPositions(atomCount: atomicNumbers.count)
    }
    
    // Create the molecule.
    var moleculeDesc = xTB_MoleculeDescriptor()
    moleculeDesc.atomicNumbers = atomicNumbers
    moleculeDesc.environment = environment
    moleculeDesc.netCharge = descriptor.netCharge
    moleculeDesc.netSpin = descriptor.netSpin
    moleculeDesc.positions = state.positions
    molecule = xTB_Molecule(descriptor: moleculeDesc)
    
    // Create the orbitals.
    orbitals = xTB_Orbitals()
    
    // Assign a self-reference to other objects, facilitating the hierarchical
    // public API.
    molecule.calculator = self
    orbitals.calculator = self
    loadHamiltonian(descriptor.hamiltonian)
  }
  
  /// Delete calculator object.
  deinit {
    xtb_delCalculator(&pointer)
  }
  
  /// Run a self-consistent field calculation and update the observables.
  public func singlepoint() {
    flushUpdateRecord()
    // ...
  }
}
