//
//  xTB_Calculator.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  var molecule: xTB_Molecule?
}

/// Singlepoint calculator.
public class xTB_Calculator {
  // Keep a reference to the molecule, so it is never deallocated while the
  // calculator is still in use.
  var molecule: xTB_Molecule
  
  var calc: xtb_TCalculator?
  
  /// Create new calculator object.
  public init(descriptor: xTB_CalculatorDescriptor) {
    guard let molecule = descriptor.molecule else {
      fatalError("Descriptor was incomplete.")
    }
    self.molecule = molecule
    
    calc = xtb_newCalculator()
    guard calc != nil else {
      fatalError("Could not create new xTB_Calculator.")
    }
  }
  
  /// Delete calculator object.
  deinit {
    xtb_delCalculator(&calc)
  }
  
  public enum Method {
    case forceField
    case tightBinding(Int)
  }
  
  /// Load a calculator.
  public func load(method: Method) {
    switch method {
    case .forceField:
      xtb_loadGFNFF(
        molecule.environment.env, molecule.mol, calc, nil)
    case .tightBinding(let version):
      if version == 2 {
        xtb_loadGFN2xTB(
          molecule.environment.env, molecule.mol, calc, nil)
      }
    }
  }
  
  public func loadTightBinding() {
    xtb_loadGFN2xTB(
      molecule.environment.env, molecule.mol, calc, nil)
  }
  
  /// Load GFN-FF calculator.
  public func loadForceField() {
    xtb_loadGFNFF(
      molecule.environment.env, molecule.mol, calc, nil)
  }
  
  // TODO: Find an ergonomic API for point charges. Perhaps as a separate
  // class, whose gradient can be queried because it stores a strong reference
  // to the calculator object.
  
  
}
