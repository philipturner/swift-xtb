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
  /// When using GFN-FF, the positions are needed to initialize force field
  /// parameters. When using tight binding, positions can be specified
  /// after initialization.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

/// Singlepoint calculator.
public class xTB_Calculator {
  var _calculator: xtb_TCalculator
  var _molecule: xtb_TMolecule
  
  var state = State()
  var updateRecord = UpdateRecord()
  var results: xTB_Results!
  
  public init(descriptor: xTB_CalculatorDescriptor) {
    guard let calc = xtb_newCalculator() else {
      fatalError("Could not create new xTB_Calculator.")
    }
    let molecule = xTB_Molecule(descriptor: descriptor)
    _calculator = calc
    _molecule = xTB_Molecule.createObject(molecule)
    
    switch descriptor.hamiltonian {
    case .forceField:
      xtb_loadGFNFF(
        xTB_Environment._environment, _molecule, _calculator, nil)
    case .tightBinding:
      xtb_loadGFN2xTB(
        xTB_Environment._environment, _molecule, _calculator, nil)
    }
    
    let externalCharges = xTB_ExternalCharges()
    let orbitals = xTB_Orbitals(descriptor: descriptor)
    state.externalCharges = externalCharges
    state.molecule = molecule
    state.orbitals = orbitals
    
    state.externalCharges.calculator = self
    state.molecule.calculator = self
    state.orbitals.calculator = self
  }
  
  deinit {
    xtb_delMolecule(&_molecule)
    xtb_delCalculator(&_calculator)
  }
}
