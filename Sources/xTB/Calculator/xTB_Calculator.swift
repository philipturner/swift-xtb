//
//  xTB_Calculator.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

public enum xTB_Hamiltonian {
  /// GFN-FF
  case forceField
  
  /// GFN2-xTB
  case tightBinding
}

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. The backing calculation environment.
  public var environment: xTB_Environment?
  
  /// Required. The parametrized method for evaluating forces.
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
  var environment: xTB_Environment
  var _calculator: xtb_TCalculator
  var _molecule: xtb_TMolecule
  
  var state = State()
  var updateRecord = UpdateRecord()
  var results: xTB_Results?
  
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
    _calculator = calc
    
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
    moleculeDesc.positions = descriptor.positions
    molecule = xTB_Molecule(descriptor: moleculeDesc)
    
    // Create the orbitals.
    orbitals = xTB_Orbitals()
    
    // Assign a self-reference to other objects, facilitating the hierarchical
    // public API.
    molecule.calculator = self
    orbitals.calculator = self
    loadHamiltonian(descriptor.hamiltonian)
  }
  
  deinit {
    xtb_delMolecule(&_molecule)
    xtb_delCalculator(&_calculator)
  }
  
  func loadHamiltonian(_ hamiltonian: xTB_Hamiltonian) {
    switch hamiltonian {
    case .forceField:
      xtb_loadGFNFF(
        environment.pointer,
        molecule.pointer,
        self.pointer,
        nil)
    case .tightBinding:
      xtb_loadGFN2_xTB(
        environment.pointer,
        molecule.pointer,
        self.pointer,
        nil)
    }
  }
  
  /// Run a self-consistent field calculation and update the observables.
  public func singlepoint() {
    flushUpdateRecord()
    // ...
  }
}
